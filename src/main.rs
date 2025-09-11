use two_round_fully_adaptive::helpers::rand_scalar;
use two_round_fully_adaptive::polynomial::*;
use two_round_fully_adaptive::sig_setup::SigParams;

fn main() {
    let n = 128;
    let d = 1 << 10;
    let sp = SigParams::setup(n, d, "MySigScheme");
    // f(X) = a0 + a1 X + a2 X^2
    let poly = Polynomial::random(100, &sp.pp);
    let z = rand_scalar();

    // Evaluate and prove
    let (y, proof, public, _) = poly.prove_eval(z, &sp.pp);

    // Verify
    let ok = verify_eval(&public, &proof, &sp.pp);
    println!(
        "Polynomial degree {} evaluated at z: verify = {ok}",
        poly.degree()
    );
}
