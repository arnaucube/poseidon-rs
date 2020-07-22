use blake2::digest::{Input, VariableOutput};
use blake2::VarBlake2b;

extern crate num;
extern crate num_bigint;
extern crate num_traits;
use num_bigint::{BigInt, Sign};
use num_traits::{One, Zero};

const SEED: &str = "poseidon";
const NROUNDSF: usize = 8;
const NROUNDSP: usize = 57;
const T: usize = 6;

#[derive(Debug)]
pub struct Constants {
    c: Vec<BigInt>,
    m: Vec<Vec<BigInt>>,
}

pub fn generate_constants() -> Constants {
    let r: BigInt = BigInt::parse_bytes(
        b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10,
    )
    .unwrap();

    let c = get_pseudo_random(&r, format!("{}{}", SEED, "_constants"), NROUNDSF + NROUNDSP);
    let m = get_mds(&r);
    Constants { c: c, m: m }
}

pub fn get_pseudo_random(r: &BigInt, seed: String, n: usize) -> Vec<BigInt> {
    let mut hasher = VarBlake2b::new(32).unwrap();
    hasher.input(seed.as_bytes());
    let mut h = hasher.vec_result();

    let mut res: Vec<BigInt> = Vec::new();
    while res.len() < n {
        let new_n: BigInt = modulus(&BigInt::from_bytes_le(Sign::Plus, &h), &r);
        res.push(new_n);
        let mut hasher = VarBlake2b::new(32).unwrap();
        hasher.input(h);
        h = hasher.vec_result();
    }
    res
}

pub fn nonce_to_string(n: usize) -> String {
    let mut r = format!("{}", n);
    while r.len() < 4 {
        r = format!("0{}", r);
    }
    r
}

pub fn get_mds(r: &BigInt) -> Vec<Vec<BigInt>> {
    let mut nonce = 0;
    let mut cauchy_matrix = get_pseudo_random(
        r,
        format!("{}_matrix_{}", SEED, nonce_to_string(nonce)),
        T * 2,
    );
    while !check_all_different(&cauchy_matrix) {
        nonce = nonce + 1;
        cauchy_matrix = get_pseudo_random(
            r,
            format!("{}_matrix_{}", SEED, nonce_to_string(nonce)),
            T * 2,
        );
    }
    let mut m: Vec<Vec<BigInt>> = Vec::new();
    for i in 0..T {
        let mut mi: Vec<BigInt> = Vec::new();
        for j in 0..T {
            mi.push(modinv(
                &modulus(&(&cauchy_matrix[i] - &cauchy_matrix[T + j]), &r),
                &r,
            ));
        }
        m.push(mi);
    }
    m
}

pub fn check_all_different(v: &Vec<BigInt>) -> bool {
    let zero: BigInt = Zero::zero();
    for i in 0..v.len() {
        if v[i] == zero {
            return false;
        }
        for j in i + 1..v.len() {
            if v[i] == v[j] {
                return false;
            }
        }
    }
    true
}

pub fn modulus(a: &BigInt, m: &BigInt) -> BigInt {
    ((a % m) + m) % m
}

pub fn modinv(a: &BigInt, q: &BigInt) -> BigInt {
    let mut mn = (q.clone(), a.clone());
    let mut xy: (BigInt, BigInt) = (Zero::zero(), One::one());

    let big_zero: BigInt = Zero::zero();
    while mn.1 != big_zero {
        xy = (xy.1.clone(), xy.0 - (mn.0.clone() / mn.1.clone()) * xy.1);
        mn = (mn.1.clone(), modulus(&mn.0, &mn.1));
    }

    while xy.0 < Zero::zero() {
        xy.0 = modulus(&xy.0, q);
    }
    xy.0
}

fn main() {
    let c = generate_constants();
    println!("let c_str: Vec<&str> = vec![");
    for i in 0..c.c.len() {
        println!("  {:?},", c.c[i].to_string());
    }
    println!("];\n");
    println!("let m_str: Vec<Vec<&str>> = vec![");
    for i in 0..c.m.len() {
        println!("  vec![");
        for j in 0..c.m[i].len() {
            println!("      {:?},", c.m[i][j].to_string());
        }
        println!("  ],");
    }
    println!("];\n");
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustc_hex::ToHex;

    #[test]
    fn test_blake2_version() {
        let mut hasher = VarBlake2b::new(32).unwrap();
        hasher.input(b"poseidon_constants");
        let h = hasher.vec_result();
        assert_eq!(
            h.to_hex(),
            "e57ba154fb2c47811dc1a2369b27e25a44915b4e4ece4eb8ec74850cb78e01b1"
        );
    }
}
