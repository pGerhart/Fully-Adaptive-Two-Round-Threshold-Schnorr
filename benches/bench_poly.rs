#![allow(non_snake_case)]
use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use two_round_fully_adaptive::helpers::{next_pow2, rand_scalar};
use two_round_fully_adaptive::polynomial::{Params, Polynomial, verify_eval};

fn bench_poly_eval(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_eval");

    // degrees = 2^4, 2^6, ..., 2^16
    for exp in (4..=16).step_by(1) {
        let degree = 1usize << exp;
        let n = degree + 1;
        let m = next_pow2(n);

        // Prepare inputs outside the timed region
        let pp = Params::setup(m, "PolyBench");
        let poly = Polynomial::random(degree, &pp);
        let z = rand_scalar();

        // --- Prove ---
        group.bench_with_input(BenchmarkId::new("prove", degree), &degree, |b, &_deg| {
            b.iter(|| {
                // pass &z and destructure 4-tuple, ignore rho_y
                let (_y, proof, public, _rho_y) = poly.prove_eval(&black_box(z), &pp);
                black_box(proof);
                black_box(public);
            });
        });

        // Precompute a proof once for verify benchmark
        let (_y_pre, proof_pre, public_pre, _rho_y_pre) = poly.prove_eval(&z, &pp);

        // --- Verify ---
        group.bench_with_input(BenchmarkId::new("verify", degree), &degree, |b, &_deg| {
            b.iter(|| {
                let ok = verify_eval(black_box(&public_pre), black_box(&proof_pre), &pp);
                black_box(ok);
            });
        });
    }

    group.finish();
}

fn bench_eval_vs_ro_eval(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_eval_functions");

    // degrees = 2^4, 2^6, ..., 2^16
    for exp in (4..=16).step_by(1) {
        let degree = 1usize << exp;
        let n = degree + 1;
        let m = next_pow2(n);

        // Setup once
        let pp = Params::setup(m, "PolyEvalFns");
        let poly = Polynomial::random(degree, &pp);
        let z = rand_scalar();
        let key = rand_scalar(); // key for on-the-fly RO-based coefficients

        // --- Stored-coeff evaluation: poly.eval(z) ---
        group.bench_with_input(
            BenchmarkId::new("eval_stored", degree),
            &degree,
            |b, &_deg| {
                b.iter(|| {
                    let y = poly.eval(&black_box(z));
                    black_box(y);
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("eval_ro_cha_cha", degree),
            &degree,
            |b, &_deg| {
                b.iter(|| {
                    let y = Polynomial::eval_with_prf_chacha12(&black_box(z), degree, &key);
                    black_box(y);
                });
            },
        );

        // --- On-the-fly RO evaluation: Polynomial::eval_with_ro(z, degree, key) ---
        group.bench_with_input(BenchmarkId::new("eval_ro", degree), &degree, |b, &_deg| {
            b.iter(|| {
                let y = Polynomial::eval_with_ro(&black_box(z), degree, &key);
                black_box(y);
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_eval_vs_ro_eval);
criterion_main!(benches);
