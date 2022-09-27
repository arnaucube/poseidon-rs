use criterion::{criterion_group, criterion_main, Criterion};

extern crate ff;
use ff::*;

use poseidon_rs::{Fr, Poseidon};

fn criterion_benchmark(c: &mut Criterion) {
    let b1: Fr = Fr::from_str(
        "12242166908188651009877250812424843524687801523336557272219921456462821518061",
    )
    .unwrap();
    let b2: Fr = Fr::from_str(
        "12242166908188651009877250812424843524687801523336557272219921456462821518061",
    )
    .unwrap();
    let mut big_arr: Vec<Fr> = Vec::new();
    big_arr.push(b1.clone());
    big_arr.push(b2.clone());
    let poseidon = Poseidon::new();

    c.bench_function("hash", |b| {
        b.iter(|| poseidon.hash(big_arr.clone()).unwrap())
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
