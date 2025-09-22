#![allow(non_snake_case)]
use merlin::Transcript;
use rand::rngs::OsRng;

use two_round_fully_adaptive::helpers::{h2p, next_pow2, rand_scalar};
use two_round_fully_adaptive::rand_eval::{rel_eval_prove, rel_eval_verify};
use two_round_fully_adaptive::sig_keygen::keygen;
use two_round_fully_adaptive::sig_setup::SigParams;

fn main() {
    // pick a degree, e.g., 2^8
    let degree = 1usize << 8;
    let n = next_pow2(degree + 1);
    let domain = format!("RelEvalSingle/deg={degree}");

    // setup
    let sp = SigParams::setup(n, degree, &domain);
    let (sk, pk) = keygen(&sp);
    let C = pk.pk;
    let Cf = pk.cf;

    let h1 = h2p(&domain, "h1", 0);
    let h2 = h2p(&domain, "h2", 0);
    let z = rand_scalar();

    // --- PROVE ---
    let mut rng = OsRng;
    let (rc, poly_pub) = rel_eval_prove(
        Transcript::new(b"RelEvalSigma"),
        &mut rng,
        &sk,
        &sp,
        &C,
        &h1,
        &h2,
        &z,
    );

    // --- VERIFY ---
    let ok = rel_eval_verify(
        Transcript::new(b"RelEvalSigma"),
        &sp,
        &C,
        &Cf,
        &h1,
        &h2,
        &poly_pub.a,
        &z,
        &rc,
    );

    println!("Verify result: {ok}");
    assert!(ok);
}
