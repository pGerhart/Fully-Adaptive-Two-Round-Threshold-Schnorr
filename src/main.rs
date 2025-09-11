use bulletproofs::LinearProof;
use curve25519_dalek::{ristretto::RistrettoPoint, scalar::Scalar, traits::MultiscalarMul};
use merlin::Transcript;
use rand::{RngCore, rngs::OsRng};
use sha2::{Digest, Sha512};

// --- helpers ---
fn hash_to_point(domain: &str, i: u64) -> RistrettoPoint {
    let mut h = Sha512::new();
    h.update(domain.as_bytes());
    h.update(&i.to_le_bytes());
    let mut bytes = [0u8; 64];
    bytes.copy_from_slice(&h.finalize());
    RistrettoPoint::from_uniform_bytes(&bytes)
}
fn msm(points: &[RistrettoPoint], scalars: &[Scalar]) -> RistrettoPoint {
    RistrettoPoint::multiscalar_mul(scalars.iter(), points.iter())
}
fn rand_scalar() -> Scalar {
    let mut w = [0u8; 64];
    OsRng.fill_bytes(&mut w);
    Scalar::from_bytes_mod_order_wide(&w)
}

// --- params ---
struct Params {
    g0: RistrettoPoint,     // F (encodes y)
    g: Vec<RistrettoPoint>, // G (encodes x)
    g_rho: RistrettoPoint,  // B (blinding base)
}
fn setup(n: usize, domain: &str) -> Params {
    let g0 = hash_to_point(&format!("{domain}/g0"), 0);
    let g_rho = hash_to_point(&format!("{domain}/g_rho"), 0);
    let g = (0..n as u64)
        .map(|i| hash_to_point(&format!("{domain}/g"), i))
        .collect();
    Params { g0, g, g_rho }
}

fn main() {
    // dimension (power of two preferred by LinearProof)
    let n = 8usize;
    let pp = setup(n, "LogDotPf");

    // public vector a
    let a: Vec<Scalar> = (0..n).map(|_| rand_scalar()).collect();

    // witness (x, rho), with y = <x, a>
    let x: Vec<Scalar> = (0..n).map(|_| rand_scalar()).collect();
    let rho = rand_scalar();
    let y = x
        .iter()
        .zip(a.iter())
        .fold(Scalar::ZERO, |acc, (xi, ai)| acc + xi * ai);

    // P = <x,G> + y*g0 + rho*g_rho
    let P = msm(&pp.g, &x) + pp.g0 * y + pp.g_rho * rho;

    // --- PROVE ---
    let mut t_prove = Transcript::new(b"LogDotPf");
    let mut rng = OsRng;

    // LinearProof expects:
    //   C: &CompressedRistretto  (the commitment)
    //   r: Scalar                (blinding opening for B)
    //   a_vec: Vec<Scalar>       (secret vector = x)
    //   b_vec: Vec<Scalar>       (public vector = a)
    //   G_vec: Vec<RistrettoPoint> (bases for a_vec)
    //   F: &RistrettoPoint       (base for c = y)
    //   B: &RistrettoPoint       (blinding base)
    let proof = LinearProof::create(
        &mut t_prove,
        &mut rng,
        &P.compress(),
        rho,
        x.clone(),
        a.clone(),
        pp.g.clone(),
        &pp.g0,
        &pp.g_rho,
    )
    .expect("prove");

    // --- VERIFY ---
    let mut t_verify = Transcript::new(b"LogDotPf");
    proof
        .verify(
            &mut t_verify,
            &P.compress(),
            &pp.g,
            &pp.g0,
            &pp.g_rho,
            a, // public vector
        )
        .expect("verify");

    println!("LinearProof verified ✅");
}
