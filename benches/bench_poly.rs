use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use two_round_fully_adaptive::helpers::{next_pow2, rand_scalar};
use two_round_fully_adaptive::polynomial::{Params, Polynomial};

fn bench_poly_eval(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_eval_prove");

    // degrees = 2^4, 2^6, ..., 2^16
    for exp in (4..=16).step_by(2) {
        let degree = 1 << exp; // polynomial degree d
        let n = degree + 1; // number of coeffs
        let m = next_pow2(n); // required by LinearProof

        // --- prepare inputs OUTSIDE the hot loop ---
        let poly = Polynomial::random(degree);
        let z = rand_scalar();
        let pp = Params::setup(m, "PolyBench"); // reuse the same params for this degree

        group.bench_with_input(BenchmarkId::from_parameter(degree), &degree, |b, &_deg| {
            b.iter(|| {
                // benchmark only proving
                let (_y, _proof, _public) = poly.prove_eval(black_box(z), &pp);
                // black_box(_proof); // uncomment if you want to keep the proof live
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_poly_eval);
criterion_main!(benches);
