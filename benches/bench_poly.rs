use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use curve25519_dalek::scalar::Scalar;
use two_round_fully_adaptive::helpers::rand_scalar;
use two_round_fully_adaptive::polynomial::Polynomial; // your module

fn bench_poly_eval(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_eval_prove");

    // degrees = 2^4, 2^6, ..., 2^16
    for exp in (4..=16).step_by(2) {
        let degree = 1 << exp;

        // prepare inputs outside of the benchmark loop
        let poly = Polynomial::random(degree);
        let z = rand_scalar();

        group.bench_with_input(BenchmarkId::from_parameter(degree), &degree, |b, &_deg| {
            b.iter(|| {
                // benchmark only proving
                let (_y, _proof, _public) = poly.prove_eval(black_box(z), "PolyBench");
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_poly_eval);
criterion_main!(benches);
