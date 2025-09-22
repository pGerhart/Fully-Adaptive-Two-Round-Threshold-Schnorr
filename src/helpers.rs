use curve25519_dalek::{
    ristretto::RistrettoPoint,
    scalar::Scalar,
    traits::{MultiscalarMul, VartimeMultiscalarMul},
};
use rand::RngCore;
use rand::rngs::OsRng;

use sha2::{Digest, Sha512};

/// Derive a Ristretto generator deterministically from (domain, tag, i).
#[inline]
pub fn h2p(domain: &str, tag: &str, i: u64) -> RistrettoPoint {
    let mut h = Sha512::new();
    h.update(domain.as_bytes());
    h.update(tag.as_bytes());
    h.update(&i.to_le_bytes());
    let mut out = [0u8; 64];
    out.copy_from_slice(&h.finalize());
    RistrettoPoint::from_uniform_bytes(&out)
}

// --------------------------------------------------------------------------------------
// MSM (multi-scalar multiplication)
// --------------------------------------------------------------------------------------

/// Constant-time MSM (side-channel hardening).
#[inline]
pub fn msm_ct(points: &[RistrettoPoint], scalars: &[Scalar]) -> RistrettoPoint {
    // Note: takes iterators; this version is constant-time but slower.
    RistrettoPoint::multiscalar_mul(scalars.iter(), points.iter())
}

#[inline]
pub fn msm_vt(points: &[RistrettoPoint], scalars: &[Scalar]) -> RistrettoPoint {
    RistrettoPoint::vartime_multiscalar_mul(scalars, points)
}

#[inline]
pub fn msm(scalars: &[Scalar], points: &[RistrettoPoint]) -> RistrettoPoint {
    msm_vt(points, scalars)
}

// --------------------------------------------------------------------------------------
// Small utils
// --------------------------------------------------------------------------------------

#[inline]
pub fn next_pow2(n: usize) -> usize {
    n.next_power_of_two()
}

#[inline]
pub fn rand_scalar() -> Scalar {
    let mut w = [0u8; 64];
    OsRng.fill_bytes(&mut w);
    Scalar::from_bytes_mod_order_wide(&w)
}

// --------------------------------------------------------------------------------------
// In-place padding / powers / evaluation
// --------------------------------------------------------------------------------------

/// Ensure `v.len() >= m` by padding zeros (does not shrink).
#[inline]
pub fn pad_zeros_in_place(v: &mut Vec<Scalar>, m: usize) {
    if v.len() < m {
        v.resize(m, Scalar::ZERO);
    }
}

/// Return a new padded vector.
#[inline]
pub fn pad_zeros(mut v: Vec<Scalar>, m: usize) -> Vec<Scalar> {
    pad_zeros_in_place(&mut v, m);
    v
}

/// Fill `powers` with `[1, z, z^2, ..., z^(n-1)]`, reusing its capacity.
#[inline]
pub fn vandermonde(powers: &mut Vec<Scalar>, z: Scalar, n: usize) {
    powers.clear();
    powers.reserve(n.saturating_sub(powers.capacity()));
    let mut pow = Scalar::ONE;
    for _ in 0..n {
        powers.push(pow);
        pow *= z;
    }
}

/// Compute `y = <coeffs, [1, z, z^2, ...]>`.
#[inline]
pub fn dot_with_vandermonde(coeffs: &[Scalar], z: Scalar) -> Scalar {
    let mut pow = Scalar::ONE;
    let mut y = Scalar::ZERO;
    for &c in coeffs {
        y += c * pow;
        pow *= z;
    }
    y
}

/// Horner's method: `y = ((((a_d) z + a_{d-1}) z + ...) z + a_0)`.
#[inline]
pub fn horner_eval(coeffs: &[Scalar], z: Scalar) -> Scalar {
    coeffs
        .iter()
        .rev()
        .fold(Scalar::ZERO, |acc, &c| acc * z + c)
}

/// Evaluate `y` and optionally write the Vandermonde powers in-place.
/// When `powers_out` is `Some`, we fill it with `[1, z, ..., z^(n-1)]`
/// while computing `y` in a single pass (cache-friendly).
#[inline]
pub fn eval_and_maybe_powers(
    coeffs: &[Scalar],
    z: Scalar,
    powers_out: Option<&mut Vec<Scalar>>,
) -> Scalar {
    match powers_out {
        Some(p) => {
            p.clear();
            p.reserve(coeffs.len().saturating_sub(p.capacity()));
            let mut pow = Scalar::ONE;
            let mut y = Scalar::ZERO;
            for &c in coeffs {
                p.push(pow);
                y += c * pow;
                pow *= z;
            }
            y
        }
        None => horner_eval(coeffs, z),
    }
}

/// Return a **new** powers vector and `y`. Prefer the in-place version above in hot paths.
#[inline]
pub fn powers_and_eval(coeffs: &[Scalar], z: Scalar) -> (Vec<Scalar>, Scalar) {
    let mut a = Vec::with_capacity(coeffs.len());
    let y = eval_and_maybe_powers(coeffs, z, Some(&mut a));
    (a, y)
}
