#![allow(non_snake_case)]
use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use two_round_fully_adaptive::helpers::{next_pow2, rand_scalar};
use two_round_fully_adaptive::polynomial::{Params, Polynomial, verify_eval};

fn bench_poly_eval(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_eval");

    // degrees = 2^4, 2^6, ..., 2^16
    for exp in (4..=16).step_by(2) {
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

criterion_group!(benches, bench_poly_eval);
criterion_main!(benches);
