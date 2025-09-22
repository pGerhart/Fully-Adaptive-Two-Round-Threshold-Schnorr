# RelEval Prototype

This repository contains a prototype implementation accompanying the paper  
**Fully-Adaptive Two-Round Threshold Schnorr Signatures from DDH**.


## Overview

We implement the cryptographic core of our threshold signature scheme, which combines:

- **Polynomial evaluation**: signing randomness is derived by evaluating a secret polynomial at a given point.
- **Proof of correct evaluation**: a polynomial commitment scheme ensures that the evaluation was done honestly.
- **Well-formedness of the first-round message**: proved by combining the polynomial commitment with a Chaum–Pedersen style proof.



## Cryptographic building blocks

- **Polynomial commitment scheme**  
  - Commitments: Pedersen commitments \(C_f\) to the polynomial coefficients.  
  - Proofs: logarithmic-size inner-product arguments (from the [dalek-cryptography/bulletproofs](https://github.com/dalek-cryptography/bulletproofs) crate) for proving correct evaluation.

- **Chaum–Pedersen proof**  
  For showing correctness of the first-round message \(R\) with respect to a signer’s public key  
  \[
  pk = (g^x h^w v^u, C_f)
  \]
  and the evaluation point \(z\).



## Code structure

- **`sig_setup.rs` / `sig_keygen.rs`**  
  Parameter setup and partial public key generation.

- **`rand_eval.rs`**  
  Core RelEval protocol: prover (`rel_eval_prove`) and verifier (`rel_eval_verify`).

- **`benches/bench_com.rs`**  
  Criterion benchmarks for proving and verifying correctness of the first-round message at different polynomial degrees.

- **`main.rs`**  
  A simple round-trip example: generate a proof, evaluate, and verify.



## Requirements

- Rust ≥ 1.75
- Cargo
- [criterion](https://crates.io/crates/criterion) (for benchmarks)



## Usage

Build the project:

```bash
cargo build
```

Run the round-trip demo (from main.rs):
```bash
cargo run 
```

Run benchmarks:
```bash
cargo bench 
```