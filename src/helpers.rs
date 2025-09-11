use bulletproofs::LinearProof;
use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
    traits::MultiscalarMul,
};
use merlin::Transcript;
use rand::{RngCore, rngs::OsRng};
use sha2::{Digest, Sha512};

// ---------- helpers ----------
pub fn h2p(domain: &str, tag: &str, i: u64) -> RistrettoPoint {
    let mut h = Sha512::new();
    h.update(domain.as_bytes());
    h.update(tag.as_bytes());
    h.update(&i.to_le_bytes());
    let mut b = [0u8; 64];
    b.copy_from_slice(&h.finalize());
    RistrettoPoint::from_uniform_bytes(&b)
}

pub fn msm(points: &[RistrettoPoint], scalars: &[Scalar]) -> RistrettoPoint {
    RistrettoPoint::multiscalar_mul(scalars.iter(), points.iter())
}

pub fn pad_zeros(mut v: Vec<Scalar>, m: usize) -> Vec<Scalar> {
    v.resize(m, Scalar::ZERO);
    v
}

pub fn next_pow2(n: usize) -> usize {
    n.next_power_of_two()
}

pub fn rand_scalar() -> Scalar {
    let mut w = [0u8; 64];
    OsRng.fill_bytes(&mut w);
    Scalar::from_bytes_mod_order_wide(&w)
}
pub fn vandermonde(z: Scalar, d: usize) -> Vec<Scalar> {
    let mut a = Vec::with_capacity(d + 1);
    let mut pow = Scalar::ONE;
    for _ in 0..=d {
        a.push(pow);
        pow *= z;
    }
    a
}

/// Horner's method: evaluates f(z) where coeffs = [x0, x1, ..., xd].
/// Returns y = x_d; for i=d-1..0: y = y*z + x_i.
pub fn horner_eval(coeffs: &[Scalar], z: Scalar) -> Scalar {
    coeffs
        .iter()
        .rev()
        .fold(Scalar::ZERO, |acc, &c| acc * z + c)
}

/// Compute [1, z, z^2, ..., z^d] and y = <coeffs, powers> in a single pass.
/// More cache-friendly than first building powers and then doing a dot-product.
pub fn powers_and_eval(coeffs: &[Scalar], z: Scalar) -> (Vec<Scalar>, Scalar) {
    let n = coeffs.len();
    let mut a = Vec::with_capacity(n);
    let mut pow = Scalar::ONE;
    let mut y = Scalar::ZERO;
    for &c in coeffs {
        a.push(pow); // current power
        y += c * pow; // accumulate f(z)
        pow *= z; // next power
    }
    (a, y)
}
