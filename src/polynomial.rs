#![allow(non_snake_case)]
use crate::helpers::{
    coeff_at, coeff_at_chacha12, coeff_from_state, h2p, msm, next_pow2, rand_scalar,
};
use bulletproofs::LinearProof;
use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
};
use merlin::{Transcript, TranscriptRng};
use rand::rngs::OsRng;
use rayon::prelude::*;
use sha2::{Digest, Sha512};

/// Parameters (bases) for Pedersen commitments.
#[derive(Clone, Debug)]
pub struct Params {
    pub g0: RistrettoPoint,     // base for y
    pub g_rho: RistrettoPoint,  // blinding base
    pub g: Vec<RistrettoPoint>, // bases for coeffs (len = m)
}

impl Params {
    pub fn setup(m: usize, domain: &str) -> Self {
        let m = next_pow2(m);
        let g0: RistrettoPoint = h2p(domain, "g_0", 0);
        let g_rho: RistrettoPoint = h2p(domain, "g_rho", 0);
        let g = (0..m as u64).map(|i| h2p(domain, "g", i)).collect();
        Self { g0, g_rho, g }
    }

    /// Like `setup`, but lets the caller *fix* g0. Useful to share the
    /// same `g` used by the signature scheme for `g^y` in the IPA.
    pub fn setup_with_g0(m: usize, domain: &str, g0: RistrettoPoint) -> Self {
        let m = next_pow2(m);
        let g_rho = h2p(domain, "g_rho", 0);
        let g = (0..m as u64).map(|i| h2p(domain, "g", i)).collect();
        Self { g0, g_rho, g }
    }
}

/// A univariate polynomial f(X) = a_0 + a_1 X + ... + a_d X^d.
#[derive(Clone, Debug)]
pub struct Polynomial {
    coeffs: Vec<Scalar>,    // length = d+1
    pub rho_f: Scalar,      // blinding for C_f
    pub cf: RistrettoPoint, // commitment to coeffs
    x_padded: Vec<Scalar>,  // cached padded coeffs (len = m)
    m: usize,               // cache m = pp.g.len()
}

// Build x_padded and m at keygen:
impl Polynomial {
    pub fn random(degree: usize, pp: &Params) -> Self {
        let coeffs: Vec<Scalar> = (0..=degree).map(|_| rand_scalar()).collect();
        let rho_f = rand_scalar();

        let n = coeffs.len();
        let m = pp.g.len();
        debug_assert!(m.is_power_of_two());
        debug_assert!(m >= n, "Params.g too short for polynomial degree");

        // cache padded coeffs once
        let mut x_padded = Vec::with_capacity(m);
        x_padded.extend_from_slice(&coeffs);
        if m > n {
            x_padded.resize(m, Scalar::ZERO);
        }

        // compute Cf once
        let cf = msm(&x_padded, &pp.g) + pp.g_rho * rho_f;

        Self {
            coeffs,
            rho_f,
            cf,
            x_padded,
            m,
        }
    }

    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    // stored-coeffs evaluator
    pub fn eval(&self, z: &Scalar) -> Scalar {
        let coeffs = &self.coeffs;

        const PAR_THRESHOLD: usize = 64;
        const B: usize = 32;

        if coeffs.len() < PAR_THRESHOLD {
            return coeffs
                .iter()
                .rev()
                .fold(Scalar::ZERO, |acc, &c| acc * z + c);
        }

        // Precompute z^k for k=0..=B
        let mut z_pows = [Scalar::ONE; B + 1];
        for k in 1..=B {
            z_pows[k] = z_pows[k - 1] * z;
        }

        // 1) Parallel: evaluate each chunk and keep its length
        let blocks: Vec<(Scalar, usize)> = coeffs
            .par_chunks(B)
            .map(|chunk| {
                let mut acc = Scalar::ZERO;
                for &c in chunk.iter().rev() {
                    acc = acc * z + c;
                }
                (acc, chunk.len())
            })
            .collect();

        // 2) Combine with z^{len(block)} (NOT constant z^B)
        blocks
            .into_iter()
            .rev()
            .fold(Scalar::ZERO, |acc, (bval, blen)| acc * z_pows[blen] + bval)
    }

    // RO-on-the-fly evaluator
    pub fn eval_with_ro(z: &Scalar, degree: usize, key: &Scalar) -> Scalar {
        const PAR_THRESHOLD: usize = 64;
        const B: usize = 32;

        // Build the pre-absorbed SHA-512 state once: DST || key
        const DST: &[u8] = b"poly-coeff";
        let mut base = Sha512::new();
        base.update(DST);
        base.update(key.to_bytes());

        // Small degrees: sequential Horner using the cloned state trick.
        if degree < PAR_THRESHOLD {
            return (0..=degree)
                .rev()
                .fold(Scalar::ZERO, |acc, i| acc * z + coeff_from_state(&base, i));
        }

        // Precompute z^k for k=0..=B (for tail block combining)
        let mut z_pows = [Scalar::ONE; B + 1];
        for k in 1..=B {
            z_pows[k] = z_pows[k - 1] * z;
        }

        // Ordered block starts for stable combine
        let starts: Vec<usize> = (0..=degree).step_by(B).collect();

        let blocks: Vec<(Scalar, usize)> = starts
            .into_par_iter()
            .map(|start| {
                let end = core::cmp::min(start + B, degree + 1);
                let blen = end - start;

                // Derive block coefficients once (forward), then reversed Horner.
                // Keep on stack for small B.
                let mut buf: [Scalar; B] = [Scalar::ZERO; B];
                for (off, i) in (start..end).enumerate() {
                    buf[off] = coeff_from_state(&base, i);
                }

                // Fast path: full block unrolled-ish
                let mut acc = Scalar::ZERO;
                for j in (0..blen).rev() {
                    acc = acc * z + unsafe { *buf.get_unchecked(j) };
                }
                (acc, blen)
            })
            .collect();

        // Combine with z^{len(block)} (correct for tail)
        blocks
            .into_iter()
            .rev()
            .fold(Scalar::ZERO, |acc, (bval, blen)| acc * z_pows[blen] + bval)
    }

    pub fn eval_with_prf_chacha12(z: &Scalar, degree: usize, key: &Scalar) -> Scalar {
        const PAR_THRESHOLD: usize = 64;
        const B: usize = 32;

        // Small degrees: sequential Horner
        if degree < PAR_THRESHOLD {
            return (0..=degree)
                .rev()
                .fold(Scalar::ZERO, |acc, i| acc * z + coeff_at_chacha12(key, i));
        }

        // Precompute z^k for k=0..=B (correct tail combination)
        let mut z_pows = [Scalar::ONE; B + 1];
        for k in 1..=B {
            z_pows[k] = z_pows[k - 1] * z;
        }

        // Ordered block starts for stable combine
        let starts: Vec<usize> = (0..=degree).step_by(B).collect();

        // Parallel per-block: derive coeffs, then reversed Horner in the block
        let blocks: Vec<(Scalar, usize)> = starts
            .into_par_iter()
            .map(|start| {
                let end = core::cmp::min(start + B, degree + 1);
                let blen = end - start;

                // derive block coeffs into a small stack buffer
                let mut buf: [Scalar; B] = [Scalar::ZERO; B];
                for (off, i) in (start..end).enumerate() {
                    buf[off] = coeff_at_chacha12(key, i);
                }

                let mut acc = Scalar::ZERO;
                for j in (0..blen).rev() {
                    // SAFETY: j < blen <= B
                    acc = acc * z + unsafe { *buf.get_unchecked(j) };
                }
                (acc, blen)
            })
            .collect();

        // Combine blocks with z^{len(block)}
        blocks
            .into_iter()
            .rev()
            .fold(Scalar::ZERO, |acc, (bval, blen)| acc * z_pows[blen] + bval)
    }

    /// Commit to coefficients with randomness rho_f.
    pub fn commit(&self, pp: &Params, rho_f: Scalar) -> RistrettoPoint {
        let n = self.coeffs.len();
        let m = pp.g.len();
        debug_assert!(m.is_power_of_two());
        debug_assert!(m >= n, "Params.g too short for polynomial degree");

        // Build a padded scalar vector once for MSM
        let mut scalars = Vec::with_capacity(m);
        scalars.extend_from_slice(&self.coeffs);
        if m > n {
            scalars.resize(m, Scalar::ZERO);
        }

        msm(&scalars, &pp.g) + pp.g_rho * rho_f
    }

    /// Produce (y, proof, public inputs) for f(z), reusing cached Cf and rho_f.
    pub fn prove_eval(&self, z: &Scalar, pp: &Params) -> (Scalar, LinearProof, PublicEval, Scalar) {
        // Safety checks
        let n = self.coeffs.len();
        let m = pp.g.len();
        debug_assert!(m == self.m, "pp.m changed since keygen");
        debug_assert!(m.is_power_of_two() && m >= next_pow2(n));

        // Build powers a(z) and y in one pass
        let mut a = Vec::with_capacity(m);
        let mut pow = Scalar::ONE;
        let mut y = Scalar::ZERO;
        for &c in &self.coeffs {
            a.push(pow);
            y += c * pow;
            pow *= z;
        }
        if m > n {
            a.resize(m, Scalar::ZERO);
        }

        // Make transcript; derive rho_y from it (faster, and FS-tied)
        let mut t = Transcript::new(b"PolyEval");
        t.append_u64(b"m", m as u64);
        t.append_message(b"g0", pp.g0.compress().as_bytes());
        t.append_message(b"g_rho", pp.g_rho.compress().as_bytes());
        t.append_message(b"z", z.as_bytes());
        let mut tr_rng: TranscriptRng = t.build_rng().finalize(&mut OsRng);

        let rho_y = Scalar::random(&mut tr_rng);
        let Cy = pp.g0 * y + pp.g_rho * rho_y;

        // Combine with cached Cf
        let Cf = self.cf;
        let P = (Cf + Cy).compress();
        let rho = self.rho_f + rho_y;

        // Create proof (must pass owned Vecs; clone cached x_padded)
        let proof = LinearProof::create(
            &mut t,
            &mut tr_rng,
            &P,
            rho,
            self.x_padded.clone(), // cheap clone vs rebuild/resize each time
            a.clone(),             // pass a to the proof; keep original for PublicEval
            pp.g.clone(),
            &pp.g0,
            &pp.g_rho,
        )
        .expect("LinearProof::create");

        let public = PublicEval {
            Cf,
            Cy,
            z: z.clone(),
            a,
            P,
        };
        (y, proof, public, rho_y)
    }
}

/// Public data for polynomial evaluation proofs.
pub struct PublicEval {
    pub Cf: RistrettoPoint,
    pub Cy: RistrettoPoint,
    pub z: Scalar,
    pub a: Vec<Scalar>,         // padded [1,z,...,z^d,0,...]
    pub P: CompressedRistretto, // Cf + Cy
}

/// Verify a polynomial evaluation proof.
pub fn verify_eval(public: &PublicEval, proof: &LinearProof, pp: &Params) -> bool {
    let mut tr = Transcript::new(b"PolyEval");

    // must match prover's absorbs EXACTLY and in the same order
    let m = public.a.len(); // safer than pp.g.len() in case of mismatch
    tr.append_u64(b"m", m as u64);
    tr.append_message(b"g0", pp.g0.compress().as_bytes());
    tr.append_message(b"g_rho", pp.g_rho.compress().as_bytes());
    tr.append_message(b"z", public.z.as_bytes());

    proof
        .verify(
            &mut tr,
            &public.P,
            &pp.g,
            &pp.g0,
            &pp.g_rho,
            public.a.clone(),
        )
        .is_ok()
}
#[cfg(test)]
mod tests {
    use super::*;
    use rand::RngCore;
    use rayon::iter::ParallelBridge; // for par_bridge()
    use rayon::prelude::*;

    /// Helper: make a Polynomial whose coeffs are derived from the RO:
    /// a_i = coeff_at(key, i), for i = 0..=degree.
    fn poly_from_ro_coeffs(degree: usize, key: &Scalar, pp: &Params) -> Polynomial {
        // Build coeffs deterministically from the RO
        let coeffs: Vec<Scalar> = (0..=degree).map(|i| coeff_at(key, i)).collect();

        // Random blinding and commitment (to satisfy struct invariants)
        let rho_f = rand_scalar();

        // Cache padded coeffs and commitment, mirroring Polynomial::random
        let n = coeffs.len();
        let m = pp.g.len();
        assert!(m.is_power_of_two());
        assert!(m >= n, "Params.g too short for polynomial degree");

        let mut x_padded = Vec::with_capacity(m);
        x_padded.extend_from_slice(&coeffs);
        if m > n {
            x_padded.resize(m, Scalar::ZERO);
        }
        let cf = msm(&x_padded, &pp.g) + pp.g_rho * rho_f;

        Polynomial {
            coeffs,
            rho_f,
            cf,
            x_padded,
            m,
        }
    }

    /// Random Scalar using OsRng.
    fn rand_scalar_os() -> Scalar {
        rand_scalar()
    }

    #[test]
    fn eval_matches_eval_with_ro_randomized() {
        // Run multiple randomized trials; cover both small and large degrees.
        // Include some hand-picked degrees around thresholds.
        let mut degrees = vec![
            0usize, 1, 2, 7, 31, 32, 33, 63, 64, 65, 128, 255, 256, 511, 512, 1023,
        ];
        // Add a few random degrees in [0, 1200]
        let mut rng = rand::thread_rng();
        for _ in 0..8 {
            degrees.push((rng.next_u64() as usize) % 1201);
        }

        for degree in degrees {
            // Setup params with enough bases (next power-of-two >= degree+1)
            let pp = Params::setup(degree + 1, "TestDomain");

            // Random key and evaluation point
            let key = rand_scalar_os();
            let z = rand_scalar_os();

            // Build polynomial with RO-derived coefficients
            let poly = poly_from_ro_coeffs(degree, &key, &pp);

            // Evaluate both ways
            let y_coeffs = poly.eval(&z);
            let y_ro = Polynomial::eval_with_ro(&z, degree, &key);

            assert_eq!(
                y_coeffs,
                y_ro,
                "mismatch for degree={} (z={}, key={})",
                degree,
                hex::encode(z.to_bytes()),
                hex::encode(key.to_bytes())
            );
        }
    }

    #[test]
    fn eval_matches_eval_with_ro_many_random_trials() {
        // Lighter degrees, more trials for sanity.
        let trials = 50;
        let mut rng = rand::thread_rng();

        for _ in 0..trials {
            let degree = (rng.next_u64() as usize) % 300; // 0..=299
            let pp = Params::setup(degree + 1, "TestDomain");

            let key = rand_scalar();
            let z = rand_scalar();

            let poly = poly_from_ro_coeffs(degree, &key, &pp);

            let y_coeffs = poly.eval(&z);
            let y_ro = Polynomial::eval_with_ro(&z, degree, &key);

            assert_eq!(y_coeffs, y_ro, "mismatch in randomized trial");
        }
    }
}
