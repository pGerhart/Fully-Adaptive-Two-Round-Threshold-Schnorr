use bulletproofs::LinearProof;
use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
    traits::VartimeMultiscalarMul,
};
use merlin::{Transcript, TranscriptRng};
use rand::{CryptoRng, RngCore, rngs::OsRng};
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
    /// Public vector a = [1, z, z^2, ..., z^d] padded to m.
    pub a: Vec<Scalar>,
}

/// Witness for the relation.
pub struct Witness {
    /// Polynomial coefficients x_0..x_d (padded to m).
    pub x: Vec<Scalar>,
    /// Evaluation y = <x, a> = f(z).
    pub y: Scalar,
    /// Opening for P with respect to g_rho: rho = rho_f + rho_y.
    pub rho: Scalar,
}

/// A Bulletproofs linear proof proving <x, a> = y while bound to P.
pub type Proof = LinearProof;

/// Return m = next power of two >= n.
#[inline]
fn next_pow2(n: usize) -> usize {
    n.next_power_of_two()
}

/// Compute a = [1, z, z^2, ..., z^(n-1)] and y = <x_raw, a>, then pad a and x to length m.
/// This is cache-friendly and avoids a separate dot-product pass.
#[inline]
fn powers_and_eval_and_pad(
    x_raw: Vec<Scalar>,
    z: Scalar,
    m: usize,
) -> (Vec<Scalar>, Vec<Scalar>, Scalar) {
    let n = x_raw.len();
    debug_assert!(m >= n && m.is_power_of_two());

    let mut a = Vec::with_capacity(m);
    let mut y = Scalar::ZERO;
    let mut pow = Scalar::ONE;

    // First n entries: build a and accumulate y
    for xi in &x_raw {
        a.push(pow);
        y += *xi * pow;
        pow *= z;
    }
    // Pad a to m
    if m > n {
        a.resize(m, Scalar::ZERO);
    }

    // Pad x to m (reuse the owned x_raw vector)
    let mut x = x_raw;
    if m > n {
        x.resize(m, Scalar::ZERO);
    }

    (a, x, y)
}

/// Commit to the coefficient vector: C_f = <x, G> + rho_f * g_rho.
#[inline]
pub fn commit_coeffs(pp: &Params, x: &[Scalar], rho_f: Scalar) -> RistrettoPoint {
    // Constant-time MSM by default; switch to vartime if you choose in your helpers.
    RistrettoPoint::vartime_multiscalar_mul(x.iter(), pp.g.iter()) + pp.g_rho * rho_f
}

/// Commit to y: C_y = y * g0 + rho_y * g_rho.
#[inline]
pub fn commit_y(pp: &Params, y: Scalar, rho_y: Scalar) -> RistrettoPoint {
    pp.g0 * y + pp.g_rho * rho_y
}

/// Build the public and witness objects for a polynomial opening at point z.
/// Assumes `pp.g.len() == m = next_pow2(n)`. We **do not** re-derive generators here.
pub fn make_instance(
    pp: &Params,
    x_raw: Vec<Scalar>, // length n = d+1 (un-padded)
    z: Scalar,
    rho_f: Scalar,
    rho_y: Scalar,
) -> (Public, Witness) {
    let n = x_raw.len();
    let m = next_pow2(n);
    assert!(
        pp.g.len() == m,
        "Params.g must have length m = next_pow2(n); pass correctly-sized Params"
    );

    // Build powers & y and pad both vectors in one pass
    let (a, x, y) = powers_and_eval_and_pad(x_raw, z, m);

    // commitments with padded vectors
    let c_f = commit_coeffs(pp, &x, rho_f);
    let c_y = commit_y(pp, y, rho_y);
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

/// Create a proof for the instance using Bulletproofs LinearProof (moves witness to avoid cloning).
#[inline]
pub fn prove_move<R: RngCore + CryptoRng>(
    mut tr: Transcript,
    public: &Public,
    mut wit: Witness,
    rng: &mut R,
) -> Proof {
    // Bind minimal public context
    tr.append_message(b"P", public.P.as_bytes());
    tr.append_u64(b"n", public.a.len() as u64);
    tr.append_message(b"g0", public.pp.g0.compress().as_bytes());
    tr.append_message(b"g_rho", public.pp.g_rho.compress().as_bytes());

    // Merlin RNG bound to transcript + caller RNG
    let mut tr_rng: TranscriptRng = tr.build_rng().finalize(rng);

    // Move x out to avoid a clone; a must be owned Vec for API, so we clone once.
    let x_owned = std::mem::take(&mut wit.x);
    let a_owned = public.a.clone();
    let g_owned = public.pp.g.clone();

    LinearProof::create(
        &mut tr,
        &mut tr_rng,
        &public.P,
        wit.rho,
        x_owned,
        a_owned,
        g_owned,
        &public.pp.g0,
        &public.pp.g_rho,
    )
    .expect("LinearProof::create")
}
