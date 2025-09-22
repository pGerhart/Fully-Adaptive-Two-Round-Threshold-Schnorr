use crate::helpers::{h2p, msm, next_pow2, rand_scalar};
use bulletproofs::LinearProof;
use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
};
use merlin::{Transcript, TranscriptRng};
use rand::rngs::OsRng;
use rayon::prelude::*;

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

    /// Evaluate f(z) with a parallel block-Horner for large polynomials.
    pub fn eval(&self, z: &Scalar) -> Scalar {
        let coeffs = &self.coeffs;

        // Threshold where threads start to pay off (tune for your machine).
        const PAR_THRESHOLD: usize = 64;
        // Block size for baby-step/giant-step Horner (64/128 are good starting points).
        const B: usize = 32;

        if coeffs.len() < PAR_THRESHOLD {
            // Original sequential Horner (fast for small n).
            return coeffs
                .iter()
                .rev()
                .fold(Scalar::ZERO, |acc, &c| acc * z + c);
        }

        // Precompute z^B
        let zB = (0..B).fold(Scalar::ONE, |acc, _| acc * z);

        // 1) Parallel: Horner-evaluate each block of size B with base z
        let blocks: Vec<Scalar> = coeffs
            .par_chunks(B)
            .map(|chunk| chunk.iter().rev().fold(Scalar::ZERO, |acc, &c| acc * z + c))
            .collect();

        // 2) Combine block results with base z^B:
        //    ((((b_k) * z^B + b_{k-1}) * z^B + ...) * z^B + b_0)
        blocks
            .into_iter()
            .rev()
            .fold(Scalar::ZERO, |acc, b| acc * zB + b)
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
