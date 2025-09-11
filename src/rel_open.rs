use bulletproofs::LinearProof;
use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
    traits::MultiscalarMul,
};
use merlin::Transcript;
use rand::{RngCore, rngs::OsRng};
use sha2::{Digest, Sha512};

/// Public parameters (bases) for the relation.
#[derive(Clone)]
pub struct Params {
    /// Base for the scalar y (the evaluated polynomial f(z)).
    pub g0: RistrettoPoint,
    /// Base for blinding randomness (used by both C_f and C_y).
    pub g_rho: RistrettoPoint,
    /// Vector bases for polynomial coefficients x_0..x_d.
    pub g: Vec<RistrettoPoint>,
}

impl Params {
    /// Derive bases deterministically for length n = d+1, domain-separated.
    pub fn setup(n: usize, domain: &str) -> Self {
        fn h2p(domain: &str, tag: &str, i: u64) -> RistrettoPoint {
            let mut h = Sha512::new();
            h.update(domain.as_bytes());
            h.update(tag.as_bytes());
            h.update(&i.to_le_bytes());
            let mut bytes = [0u8; 64];
            bytes.copy_from_slice(&h.finalize());
            RistrettoPoint::from_uniform_bytes(&bytes)
        }
        let g0 = h2p(domain, "g0", 0);
        let g_rho = h2p(domain, "g_rho", 0);
        let g = (0..n as u64).map(|i| h2p(domain, "g", i)).collect();
        Params { g0, g_rho, g }
    }
}

/// Public inputs for the LogDotPf instance used for polynomial opening.
#[derive(Clone)]
pub struct Public {
    /// Commitment P = C_f + C_y = <x,G> + y*g0 + (rho_f+rho_y)*g_rho.
    pub P: CompressedRistretto,
    /// Bases (copied so the verifier has everything it needs).
    pub pp: Params,
    /// Public vector a = [1, z, z^2, ..., z^d].
    pub a: Vec<Scalar>,
}

/// Witness for the relation.
pub struct Witness {
    /// Polynomial coefficients x_0..x_d.
    pub x: Vec<Scalar>,
    /// Evaluation y = <x, a> = f(z).
    pub y: Scalar,
    /// Opening for P with respect to g_rho: rho = rho_f + rho_y.
    pub rho: Scalar,
}

/// A Bulletproofs linear proof proving <x, a> = y while bound to P.
pub type Proof = LinearProof;

/// Convenience: multiscalar multiplication helper.
fn msm(points: &[RistrettoPoint], scalars: &[Scalar]) -> RistrettoPoint {
    RistrettoPoint::multiscalar_mul(scalars.iter(), points.iter())
}

/// Return [1, z, z^2, ..., z^d].
pub fn vandermonde(z: Scalar, d: usize) -> Vec<Scalar> {
    let mut a = Vec::with_capacity(d + 1);
    let mut pow = Scalar::ONE;
    for _ in 0..=d {
        a.push(pow);
        pow *= z;
    }
    a
}

fn next_pow2(n: usize) -> usize {
    n.next_power_of_two()
}
fn pad_zeros(mut v: Vec<Scalar>, m: usize) -> Vec<Scalar> {
    v.resize(m, Scalar::ZERO);
    v
}

/// Sample a random Scalar.
pub fn rand_scalar() -> Scalar {
    let mut w = [0u8; 64];
    OsRng.fill_bytes(&mut w);
    Scalar::from_bytes_mod_order_wide(&w)
}

/// Commit to the coefficient vector: C_f = <x, G> + rho_f * g_rho.
pub fn commit_coeffs(pp: &Params, x: &[Scalar], rho_f: Scalar) -> RistrettoPoint {
    msm(&pp.g, x) + pp.g_rho * rho_f
}

/// Commit to y: C_y = y * g0 + rho_y * g_rho.
pub fn commit_y(pp: &Params, y: Scalar, rho_y: Scalar) -> RistrettoPoint {
    pp.g0 * y + pp.g_rho * rho_y
}

/// Build the public and witness objects for a polynomial opening at point z.
///
/// Inputs:
/// - `pp`: parameters for n=d+1
/// - `x`: coefficients (length n)
/// - `z`: challenge point
/// - `rho_f`, `rho_y`: independent blindings
///
/// Output:
/// - (`Public`, `Witness`)
pub fn make_instance(
    mut pp: Params,
    x_raw: Vec<Scalar>, // length n = d+1 (un-padded)
    z: Scalar,
    rho_f: Scalar,
    rho_y: Scalar,
) -> (Public, Witness) {
    let n = x_raw.len(); // n = d+1
    let m = next_pow2(n); // power-of-two length for LinearProof

    // Ensure we have m bases (if caller gave fewer, re-setup deterministically)
    if pp.g.len() != m {
        // keep the same domain used by caller; if you keep it external, pass it in.
        pp = Params::setup(m, "LogDotPf/poly"); // or plumb a domain param through
    }

    // Vandermonde of length n, then pad to length m
    let a_raw = vandermonde(z, n - 1);
    let a = pad_zeros(a_raw, m);

    // y = <x, a> (same if computed before padding)
    let y = x_raw
        .iter()
        .zip(a.iter())
        .fold(Scalar::ZERO, |acc, (xi, ai)| acc + xi * ai);

    // pad x to length m
    let x = pad_zeros(x_raw, m);

    // commitments with padded vectors
    let c_f = commit_coeffs(&pp, &x, rho_f); // uses m G-bases
    let c_y = commit_y(&pp, y, rho_y);
    let P = (c_f + c_y).compress();
    let rho = rho_f + rho_y;

    let public = Public {
        P,
        pp: pp.clone(),
        a,
    };
    let witness = Witness { x, y, rho };
    (public, witness)
}

/// Create a proof for the instance using Bulletproofs LinearProof.
pub fn prove(mut tr: Transcript, public: &Public, wit: &Witness) -> Proof {
    let mut rng = OsRng;
    LinearProof::create(
        &mut tr,
        &mut rng,
        &public.P,
        wit.rho,
        wit.x.clone(),
        public.a.clone(),
        public.pp.g.clone(),
        &public.pp.g0,
        &public.pp.g_rho,
    )
    .expect("LinearProof::create")
}

/// Verify a proof for the instance.
pub fn verify(mut tr: Transcript, public: &Public, proof: &Proof) -> bool {
    proof
        .verify(
            &mut tr,
            &public.P,
            &public.pp.g,
            &public.pp.g0,
            &public.pp.g_rho,
            public.a.clone(),
        )
        .is_ok()
}
