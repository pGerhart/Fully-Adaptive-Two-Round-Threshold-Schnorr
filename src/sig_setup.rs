use crate::helpers::{h2p, next_pow2};
use crate::polynomial::Params;
use curve25519_dalek::ristretto::RistrettoPoint;

/// Public parameters for the signature scheme + IPA.
/// We set t = n (your requirement).
pub struct SigParams {
    pub n: usize,          // number of signers (t = n)
    pub d: usize,          // polynomial degree bound for f
    pub g: RistrettoPoint, // scheme generator g
    pub h: RistrettoPoint, // scheme generator h
    pub v: RistrettoPoint, // scheme generator v
    pub pp: Params,        // IPA params with g0 == g
    pub domain: String,    // saved for convenience/debug
}

impl SigParams {
    /// Setup(n, t=n, d): derive (g, h, v) and IPA params with g0 = g.
    /// `domain` is used for domain separation of all bases.
    pub fn setup(n: usize, d: usize, domain: &str) -> Self {
        // signature scheme bases
        let g = h2p(domain, "sig_g", 0);
        let h = h2p(domain, "sig_h", 0);
        let v = h2p(domain, "sig_v", 0);

        // IPA length for polynomial openings: m = next_pow2(d+1)
        let m = next_pow2(d + 1);

        // make IPA params where g0 == signature g (so Cy = g^y * g_rho^rho)
        // use a sub-domain to separate the vector bases from signature bases
        let pp_domain = format!("{}/ipa", domain);
        let pp = Params::setup_with_g0(m, &pp_domain, g);

        SigParams {
            n,
            d,
            g,
            h,
            v,
            pp,
            domain: domain.to_owned(),
        }
    }
}
