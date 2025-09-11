use crate::helpers::{h2p, next_pow2, rand_scalar};
use bulletproofs::LinearProof;
use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
    traits::VartimeMultiscalarMul,
};
use merlin::{Transcript, TranscriptRng};
use rand::rngs::OsRng;

/// Parameters (bases) for Pedersen commitments.
#[derive(Clone)]
pub struct Params {
    pub g0: RistrettoPoint,     // base for y
    pub g_rho: RistrettoPoint,  // blinding base
    pub g: Vec<RistrettoPoint>, // bases for coeffs (len = m)
}
impl Params {
    pub fn setup(m: usize, domain: &str) -> Self {
        let m = next_pow2(m);
        let g0 = h2p(domain, "g0", 0);
        let g_rho = h2p(domain, "g_rho", 0);
        let g = (0..m as u64).map(|i| h2p(domain, "g", i)).collect();
        Self { g0, g_rho, g }
    }
}

/// A univariate polynomial f(X) = a_0 + a_1 X + ... + a_d X^d.
pub struct Polynomial {
    coeffs: Vec<Scalar>,   // length = d+1
    rho_f: Scalar,         // blinding for C_f
    cf: RistrettoPoint,    // commitment to coeffs
    x_padded: Vec<Scalar>, // cached padded coeffs (len = m)
    m: usize,              // cache m = pp.g.len()
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
        let cf = RistrettoPoint::vartime_multiscalar_mul(&x_padded, &pp.g) + pp.g_rho * rho_f;

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

    /// Evaluate f(z) with Horner's rule.
    pub fn eval(&self, z: Scalar) -> Scalar {
        self.coeffs
            .iter()
            .rev()
            .fold(Scalar::ZERO, |acc, &c| acc * z + c)
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

        // Fast variable-time MSM (OK for benchmarking; avoid in constant-time contexts)
        RistrettoPoint::vartime_multiscalar_mul(&scalars, &pp.g) + pp.g_rho * rho_f
    }

    /// Produce (y, proof, public inputs) for f(z), reusing cached Cf and rho_f.
    pub fn prove_eval(&self, z: Scalar, pp: &Params) -> (Scalar, LinearProof, PublicEval) {
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

        let public = PublicEval { Cf, Cy, z, a, P };
        (y, proof, public)
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
