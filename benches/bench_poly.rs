use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use two_round_fully_adaptive::helpers::{next_pow2, rand_scalar};
use two_round_fully_adaptive::polynomial::{Params, Polynomial, verify_eval};

fn bench_poly_eval(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_eval");

    for exp in (4..=16).step_by(2) {
        let degree = 1 << exp;
        let n = degree + 1;
        let m = next_pow2(n);

        let pp = Params::setup(m, "PolyBench");
        let poly = Polynomial::random(degree, &pp);
        let z = rand_scalar();

        // --- Prove ---
        group.bench_with_input(BenchmarkId::new("prove", degree), &degree, |b, &_deg| {
            b.iter(|| {
                let (_y, _proof, _public) = poly.prove_eval(black_box(z), &pp);
                // keep values alive
                black_box(&_proof);
                black_box(&_public);
            });
        });

        // Prepare a proof once for verify benchmark
        let (_y, proof, public) = poly.prove_eval(z, &pp);

        // --- Verify ---
        group.bench_with_input(BenchmarkId::new("verify", degree), &degree, |b, &_deg| {
            b.iter(|| {
                let ok = verify_eval(black_box(&public), black_box(&proof), &pp);
                black_box(ok);
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_poly_eval);
criterion_main!(benches);
