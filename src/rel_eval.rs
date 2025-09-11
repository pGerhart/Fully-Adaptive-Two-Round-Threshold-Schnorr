use crate::sig_setup::SigParams;
use curve25519_dalek::{ristretto::RistrettoPoint, scalar::Scalar, traits::VartimeMultiscalarMul};
use merlin::Transcript;
use rand::{CryptoRng, RngCore};

/// Compact Σ-proof for the relation:
///   R  = g^r * h1^w * h2^u
///   Cr = g^r * g_rho^{rho_r}
///   C  = g^x * h^w * v^u
/// Only (c, z_*) are sent; T1,T2,T3 are recomputed by the verifier.
#[derive(Clone, Debug)]
pub struct RelEvalProof {
    pub c: Scalar,
    pub z_r: Scalar,
    pub z_rho: Scalar,
    pub z_u: Scalar,
    pub z_w: Scalar,
    pub z_x: Scalar,
}

#[inline]
fn challenge_scalar(tr: &mut Transcript, label: &'static [u8]) -> Scalar {
    let mut buf = [0u8; 64];
    tr.challenge_bytes(label, &mut buf);
    Scalar::from_bytes_mod_order_wide(&buf)
}

#[inline]
fn absorb_for_challenge(
    tr: &mut Transcript,
    // commitments reconstructed (or computed by prover)
    T1: &RistrettoPoint,
    T2: &RistrettoPoint,
    T3: &RistrettoPoint,
    // public statements
    R: &RistrettoPoint,
    Cr: &RistrettoPoint,
    C: &RistrettoPoint,
) -> Scalar {
    tr.append_message(b"T1", T1.compress().as_bytes());
    tr.append_message(b"T2", T2.compress().as_bytes());
    tr.append_message(b"T3", T3.compress().as_bytes());
    tr.append_message(b"R", R.compress().as_bytes());
    tr.append_message(b"Cr", Cr.compress().as_bytes());
    tr.append_message(b"C", C.compress().as_bytes());
    challenge_scalar(tr, b"c")
}

/// Prover side (Fiat–Shamir): outputs only (c, z_*).
pub fn prove_rel_eval<RNG: RngCore + CryptoRng>(
    mut tr: Transcript,
    rng: &mut RNG,
    sp: &SigParams,
    h1: &RistrettoPoint,
    h2: &RistrettoPoint,
    C: &RistrettoPoint,
    R: &RistrettoPoint,
    Cr: &RistrettoPoint,
    r: Scalar,
    rho_r: Scalar,
    u: Scalar,
    w: Scalar,
    x: Scalar,
) -> RelEvalProof {
    // randomizers
    let k_r = Scalar::random(rng);
    let k_rho = Scalar::random(rng);
    let k_u = Scalar::random(rng);
    let k_w = Scalar::random(rng);
    let k_x = Scalar::random(rng);

    // T1 = g^{k_r} g_rho^{k_rho}
    let T1 = RistrettoPoint::vartime_multiscalar_mul([k_r, k_rho], [sp.g, sp.pp.g_rho]);

    // T2 = g^{k_r} h1^{k_w} h2^{k_u}
    let T2 = RistrettoPoint::vartime_multiscalar_mul([k_r, k_w, k_u], [sp.g, *h1, *h2]);

    // T3 = g^{k_x} h^{k_w} v^{k_u}
    let T3 = RistrettoPoint::vartime_multiscalar_mul([k_x, k_w, k_u], [sp.g, sp.h, sp.v]);

    // Fiat–Shamir challenge
    let c = absorb_for_challenge(&mut tr, &T1, &T2, &T3, R, Cr, C);

    // Responses
    let z_r = k_r + c * r;
    let z_rho = k_rho + c * rho_r;
    let z_u = k_u + c * u;
    let z_w = k_w + c * w;
    let z_x = k_x + c * x;

    RelEvalProof {
        c,
        z_r,
        z_rho,
        z_u,
        z_w,
        z_x,
    }
}

/// Verifier: recompute T’s from (c, z_*) and public data, then re-derive c’.
pub fn verify_rel_eval(
    mut tr: Transcript,
    sp: &SigParams,
    h1: &RistrettoPoint,
    h2: &RistrettoPoint,
    C: &RistrettoPoint,
    R: &RistrettoPoint,
    Cr: &RistrettoPoint,
    proof: &RelEvalProof,
) -> bool {
    let c = proof.c;

    // Recompute T1' = g^{z_r} g_rho^{z_rho} * Cr^{-c}
    let T1p = RistrettoPoint::vartime_multiscalar_mul([proof.z_r, proof.z_rho, -c], [
        sp.g,
        sp.pp.g_rho,
        *Cr,
    ]);

    // Recompute T2' = g^{z_r} h1^{z_w} h2^{z_u} * R^{-c}
    let T2p = RistrettoPoint::vartime_multiscalar_mul([proof.z_r, proof.z_w, proof.z_u, -c], [
        sp.g, *h1, *h2, *R,
    ]);

    // Recompute T3' = g^{z_x} h^{z_w} v^{z_u} * C^{-c}
    let T3p = RistrettoPoint::vartime_multiscalar_mul([proof.z_x, proof.z_w, proof.z_u, -c], [
        sp.g, sp.h, sp.v, *C,
    ]);

    // Recreate challenge and compare
    let c_prime = absorb_for_challenge(&mut tr, &T1p, &T2p, &T3p, R, Cr, C);
    c_prime == c
}
