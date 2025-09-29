#![allow(non_snake_case)]
use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use merlin::Transcript;
use rand::rngs::OsRng;

use two_round_fully_adaptive::helpers::{h2p, next_pow2, rand_scalar};
use two_round_fully_adaptive::rand_eval::{rel_eval_prove, rel_eval_verify};
use two_round_fully_adaptive::sig_keygen::keygen;
use two_round_fully_adaptive::sig_setup::SigParams;

fn bench_rel_eval(c: &mut Criterion) {
    let mut group = c.benchmark_group("rel_eval");

    // degrees = 2^4, 2^6, ..., 2^16
    for exp in (4..=16).step_by(1) {
        let degree = 1usize << exp;
        let n = next_pow2(degree + 1); // t = n in your scheme
        let domain = format!("RelEvalBench/deg={degree}");

        // --- one-time setup per degree (not in timed region) ---
        let sp = SigParams::setup(n, degree, &domain);
        let (sk, pk) = keygen(&sp);
        let C = pk.pk;
        let Cf = pk.cf;

        let h1 = h2p(&domain, "h1", 0);
        let h2 = h2p(&domain, "h2", 0);
        let z = rand_scalar();

        // Precompute a valid bundle for the verify-only benchmark
        let mut rng_pre = OsRng;
        let bundle_for_verify = rel_eval_prove(
            Transcript::new(b"RelEvalSigma"),
            &mut rng_pre,
            &sk,
            &sp,
            &C,
            &h1,
            &h2,
            &z,
        );
        // rel_eval_prove returns (RandCommitment, PublicEval)
        let (rc_pre, poly_pub_pre) = &bundle_for_verify;

        // --- PROVE benchmark ---
        group.bench_with_input(BenchmarkId::new("prove", degree), &degree, |b, &_deg| {
            let mut rng = OsRng;
            b.iter(|| {
                // measure proving end-to-end
                let bundle = rel_eval_prove(
                    Transcript::new(b"RelEvalSigma"),
                    &mut rng,
                    &sk,
                    &sp,
                    &C,
                    &h1,
                    &h2,
                    &z,
                );
                black_box(bundle);
            });
        });

        // --- VERIFY benchmark (using precomputed bundle) ---
        group.bench_with_input(BenchmarkId::new("verify", degree), &degree, |b, &_deg| {
            b.iter(|| {
                // rel_eval_verify(
                //   tr, &sp, &C, &Cf, &h1, &h2, &a: Vec<Scalar>, &z: Scalar, &RandCommitment
                // )
                let ok = rel_eval_verify(
                    Transcript::new(b"RelEvalSigma"),
                    &sp,
                    &C,
                    &Cf,
                    &h1,
                    &h2,
                    black_box(&poly_pub_pre.a),
                    black_box(&z),
                    black_box(rc_pre),
                );
                assert!(ok);
                black_box(ok);
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_rel_eval);
criterion_main!(benches);
