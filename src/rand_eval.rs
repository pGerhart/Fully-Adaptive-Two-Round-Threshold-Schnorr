// src/sig_eval.rs
use crate::helpers::msm;
use crate::polynomial::PublicEval as PolyPublic; // for constructing PublicEval on verifier
use crate::rel_eval::{RelEvalProof, prove_rel_eval, verify_rel_eval};
use crate::sig_keygen::UserSecret;
use crate::sig_setup::SigParams;
use curve25519_dalek::{ristretto::RistrettoPoint, scalar::Scalar};
use merlin::Transcript;
use rand::{CryptoRng, RngCore};

/// Compute `r = f(z)` and `R = g^r * h1^w * h2^u`.
#[inline]
pub fn compute_r(
    sk: &UserSecret,
    sp: &SigParams,
    h1: &RistrettoPoint,
    h2: &RistrettoPoint,
    z: &Scalar,
) -> (Scalar, RistrettoPoint) {
    let r = sk.poly.eval(z);
    let R = msm(&[r, sk.w, sk.u], &[sp.g, *h1, *h2]);
    (r, R)
}

/// Minimal prover bundle for rel_eval.
/// We include `a` and `P` (from the IPA public) so the verifier uses exactly what the prover used.
pub struct RandCommitment {
    pub R: RistrettoPoint,                    // session commitment
    pub Cr: RistrettoPoint,                   // evaluation commitment C_r
    pub ipa_proof: bulletproofs::LinearProof, // IPA proof for C_r
    pub sigma: RelEvalProof,                  // Σ-proof linking (R, C_r, C)
}

/// Prover: compute (r,R), prove C_r via IPA, then Σ-proof for (R,C_r,C).
/// Requires Polynomial::prove_eval_with_rho(z,&pp) -> (y, proof, public, rho_r)
pub fn rel_eval_prove<RNG: RngCore + CryptoRng>(
    tr_sigma: Transcript,
    rng: &mut RNG,
    sk: &UserSecret,
    sp: &SigParams,
    C: &RistrettoPoint,
    h1: &RistrettoPoint,
    h2: &RistrettoPoint,
    z: &Scalar,
) -> (RandCommitment, PolyPublic) {
    // 1) r and R
    let (r, R) = compute_r(sk, sp, h1, h2, z);

    // 2) IPA for C_r = g^r g_ρ^{ρ_r}
    let (r_again, ipa_proof, poly_pub, rho_r) = sk.poly.prove_eval(z, &sk.pp);
    debug_assert_eq!(r_again, r, "poly eval mismatch");

    // 3) Σ-proof ties R, C_r, C with witnesses (r, ρ_r, u, w, x)
    let sigma = prove_rel_eval(
        tr_sigma,
        rng,
        sp,
        h1,
        h2,
        C,
        &R,
        &poly_pub.Cy, // C_r
        r,
        rho_r,
        sk.u,
        sk.w,
        sk.x,
    );

    (
        RandCommitment {
            R,
            Cr: poly_pub.Cy,
            ipa_proof,
            sigma,
        },
        poly_pub,
    )
}

/// Verifier: consumes the prover’s `a` and `P` directly (no reconstruction).
/// Returns true iff IPA and Σ pass.
pub fn rel_eval_verify(
    tr_sigma: Transcript,
    sp: &SigParams,
    C: &RistrettoPoint,  // user partial public key
    Cf: &RistrettoPoint, // from UserPublic (pk.cf)
    h1: &RistrettoPoint,
    h2: &RistrettoPoint,
    a: &Vec<Scalar>,
    z: &Scalar, // public evaluation point
    com: &RandCommitment,
) -> bool {
    // Build the exact PolyPublic the prover committed to
    let poly_pub = PolyPublic {
        Cf: *Cf,
        Cy: com.Cr,
        z: *z,
        a: a.clone(),
        P: (Cf + com.Cr).compress(),
    };

    // 1) verify IPA
    let ipa_ok = crate::polynomial::verify_eval(&poly_pub, &com.ipa_proof, &sp.pp);
    if !ipa_ok {
        return false;
    }

    // 2) verify Σ proof for (R, C_r=Cy, C)
    verify_rel_eval(tr_sigma, sp, h1, h2, C, &com.R, &com.Cr, &com.sigma)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::helpers::{h2p, rand_scalar};
    use crate::sig_keygen::keygen;
    use crate::sig_setup::SigParams;
    use merlin::Transcript;
    use rand::rngs::OsRng;

    #[test]
    fn rel_eval_end_to_end_ok() {
        let n = 16usize;
        let d = 128usize;
        let domain = "RelEvalTest/v3";

        let sp = SigParams::setup(n, d, domain);
        let (sk, pk) = keygen(&sp);
        let C = pk.pk;
        let Cf = pk.cf;

        let h1 = h2p(domain, "h1", 0);
        let h2 = h2p(domain, "h2", 0);
        let z = rand_scalar();

        let mut rng = OsRng;
        let (com, polypub) = rel_eval_prove(
            Transcript::new(b"RelEvalSigma"),
            &mut rng,
            &sk,
            &sp,
            &C,
            &h1,
            &h2,
            &z,
        );

        let ok = rel_eval_verify(
            Transcript::new(b"RelEvalSigma"),
            &sp,
            &C,
            &Cf,
            &h1,
            &h2,
            &polypub.a,
            &z,
            &com,
        );
        assert!(ok, "rel_eval verification should pass");
    }

    #[test]
    fn rel_eval_sigma_tamper_fails() {
        let n = 8usize;
        let d = 64usize;
        let domain = "RelEvalTest/v3-bad";

        let sp = SigParams::setup(n, d, domain);
        let (sk, pk) = keygen(&sp);
        let C = pk.pk;
        let Cf = pk.cf;

        let h1 = h2p(domain, "h1", 1);
        let h2 = h2p(domain, "h2", 1);
        let z = rand_scalar();

        let mut rng = OsRng;
        let (mut com, polypub) = rel_eval_prove(
            Transcript::new(b"RelEvalSigma"),
            &mut rng,
            &sk,
            &sp,
            &C,
            &h1,
            &h2,
            &z,
        );

        com.sigma.c += Scalar::ONE;
        let ok = rel_eval_verify(
            Transcript::new(b"RelEvalSigma"),
            &sp,
            &C,
            &Cf,
            &h1,
            &h2,
            &polypub.a,
            &z,
            &com,
        );

        assert!(!ok, "tampered sigma proof must fail");
    }
}
