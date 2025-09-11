use merlin::Transcript;
use two_round_fully_adaptive::helpers::rand_scalar;
use two_round_fully_adaptive::polynomial::*;

fn main() {
    let pp = Params::setup(100, "PolyDomain");
    // f(X) = a0 + a1 X + a2 X^2
    let poly = Polynomial::random(100, &pp);
    let z = rand_scalar();

    // Evaluate and prove
    let (y, proof, public) = poly.prove_eval(z, &pp);

    // Verify
    let ok = verify_eval(&public, y, &proof, &pp);
    println!(
        "Polynomial degree {} evaluated at z: verify = {ok}",
        poly.degree()
    );
}
