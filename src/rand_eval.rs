// src/sig_eval.rs
use crate::polynomial::{Params as PolyParams, Polynomial, PublicEval as PolyPublic};
use crate::rel_eval::{RelEvalProof, prove_rel_eval, verify_rel_eval};
use crate::sig_keygen::UserSecret;
use crate::sig_setup::SigParams;
use curve25519_dalek::{ristretto::RistrettoPoint, scalar::Scalar, traits::VartimeMultiscalarMul};
use merlin::Transcript;
use rand::{CryptoRng, RngCore};

/// Compute `r = f(z)` and `R = g^r * h1^w * h2^u`.
///
/// Inputs:
/// - `sk`: user's long-term secret (holds w, u, x, and the polynomial f)
/// - `sp`: setup params (provides the shared base `g` and `g_ρ` via `sp.pp`)
/// - `h1`, `h2`: session-specific Ristretto points (public inputs)
/// - `z`: challenge in the scalar field
///
/// Output:
/// - `(r, R)` where r = f(z), R = g^r * h1^w * h2^u
#[inline]
pub fn compute_r(
    sk: &UserSecret,
    sp: &SigParams,
    h1: &RistrettoPoint,
    h2: &RistrettoPoint,
    z: &Scalar,
) -> (Scalar, RistrettoPoint) {
    // evaluate the nonce polynomial
    let r = sk.poly.eval(z);

    // R = g^r * h1^w * h2^u (single MSM is faster than 3 muls)
    let R = RistrettoPoint::vartime_multiscalar_mul([&r, &sk.w, &sk.u], [&sp.g, h1, h2]);
    (r, R)
}

/// Everything the prover sends for rel_eval:
/// - `poly_pub` and `ipa_proof`: IPA opening that `C_r = g^r g_ρ^{ρ_r}` with respect to `sp.pp`
/// - `sigma`: compact Σ-proof that `R` and `C` are consistent with the same witnesses
pub struct RelEvalBundle {
    pub R: RistrettoPoint,    // Public Random Commitment
    pub poly_pub: PolyPublic, // contains Cf, Cy (= C_r), z, a, P
    pub ipa_proof: bulletproofs::LinearProof,
    pub sigma: RelEvalProof, // (c, z_r, z_ρ, z_u, z_w, z_x)
}

impl RelEvalBundle {
    #[inline]
    pub fn C_r(&self) -> &RistrettoPoint {
        // Cy in your polynomial PublicEval is exactly C_r, because sp.pp.g0 == sp.g
        &self.poly_pub.Cy
    }
}

/// Prover: compute (r,R), prove C_r via IPA, then prove R is well-formed via Σ-proof.
/// - `C` is the signer's public key `g^x h^w v^u`
/// - `h1`, `h2` are session bases
/// - `z` is the evaluation point
pub fn rel_eval_prove<RNG: RngCore + CryptoRng>(
    tr_sigma: Transcript,
    rng: &mut RNG,
    sk: &UserSecret,
    sp: &SigParams,
    C: &RistrettoPoint,
    h1: &RistrettoPoint,
    h2: &RistrettoPoint,
    z: &Scalar,
) -> RelEvalBundle {
    // 1) compute r and R
    let (r, R) = compute_r(sk, sp, h1, h2, z);

    // 2) IPA: prove C_r = g^r g_ρ^{ρ_r} and bind to P = Cf + C_r (your linear proof)
    //    We need ρ_r (the blinding of C_r) for the Σ-proof.
    let (r_again, ipa_proof, poly_pub, rho_r) = sk.poly.prove_eval(*z, &sk.pp);
    debug_assert_eq!(r_again, r, "poly eval mismatch");

    // 3) Σ-proof that (R, C_r, C) are consistent with the same secrets (r, ρ_r, u, w, x)
    let sigma = prove_rel_eval(
        tr_sigma,
        rng,
        sp,
        h1,
        h2,
        C,
        &R,
        &poly_pub.Cy, // this is C_r (since g0 == g)
        r,
        rho_r,
        sk.u,
        sk.w,
        sk.x,
    );

    RelEvalBundle {
        R,
        poly_pub,
        ipa_proof,
        sigma,
    }
}

/// Verifier: check the IPA opening for `C_r` and the Σ-proof linking `R`, `C_r`, `C`.
/// Returns `true` iff both pass.
pub fn rel_eval_verify(
    tr_sigma: Transcript,
    sp: &SigParams,
    C: &RistrettoPoint,
    h1: &RistrettoPoint,
    h2: &RistrettoPoint,
    bundle: &RelEvalBundle,
) -> bool {
    // 1) verify IPA (polynomial opening)
    //    your `verify_eval` already ignores y (it recomputes via linear relation),
    //    but if it needs y, pass bundle.r. Use the same API as your module exposes.
    let ipa_ok = crate::polynomial::verify_eval(&bundle.poly_pub, &bundle.ipa_proof, &sp.pp);

    if !ipa_ok {
        return false;
    }

    // 2) verify Σ-proof
    let sigma_ok = verify_rel_eval(
        tr_sigma,
        sp,
        h1,
        h2,
        C,
        &bundle.R,
        &bundle.poly_pub.Cy, // C_r
        &bundle.sigma,
    );

    sigma_ok
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::helpers::{h2p, rand_scalar};
    use crate::sig_keygen::keygen;
    use crate::sig_setup::SigParams;
    use merlin::Transcript;
    use rand::rngs::OsRng;

    /// End-to-end: rel_eval prove + verify should succeed.
    #[test]
    fn rel_eval_end_to_end_ok() {
        // parameters: t = n (your scheme requirement), choose modest sizes for CI
        let n = 16usize;
        let d = 128usize;
        let domain = "RelEvalTest/v1";

        // global setup
        let sp = SigParams::setup(n, d, domain);

        // user keygen (includes random polynomial with cached Cf & rho_f)
        let (sk, pk) = keygen(&sp);
        let C = pk.pk; // public key point

        // session bases and evaluation point (deterministic from domain)
        let h1 = h2p(domain, "h1", 0);
        let h2 = h2p(domain, "h2", 0);
        let z = rand_scalar();

        // prover
        let mut rng = OsRng;
        let bundle = super::rel_eval_prove(
            Transcript::new(b"RelEvalSigma"),
            &mut rng,
            &sk,
            &sp,
            &C,
            &h1,
            &h2,
            &z,
        );

        // verifier
        let ok =
            super::rel_eval_verify(Transcript::new(b"RelEvalSigma"), &sp, &C, &h1, &h2, &bundle);

        assert!(ok, "rel_eval verification should pass");
    }
}
