extern crate num;
extern crate num_bigint;
extern crate num_traits;

use blake2::digest::{Input, VariableOutput};
use blake2::VarBlake2b;

use num_bigint::{BigInt, Sign};
use num_traits::{One, Zero};

const SEED: &str = "poseidon";
const NROUNDSF: usize = 8;
const NROUNDSP: usize = 57;
const T: usize = 6;

pub struct Constants {
    r: BigInt,
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
    Constants { r: r, c: c, m: m }
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

pub fn check_bigint_in_field(a: &BigInt, q: &BigInt) -> bool {
    if a >= q {
        return false;
    }
    true
}

pub fn check_bigint_array_in_field(arr: &Vec<BigInt>, q: &BigInt) -> bool {
    for a in arr {
        if !check_bigint_in_field(a, &q) {
            return false;
        }
    }
    true
}

pub struct Poseidon {
    constants: Constants,
}
impl Poseidon {
    pub fn new() -> Poseidon {
        Poseidon {
            constants: generate_constants(),
        }
    }
    pub fn ark(&self, state: &Vec<BigInt>, c: &BigInt) -> Vec<BigInt> {
        let mut new_state: Vec<BigInt> = state.clone();
        for i in 0..state.len() {
            new_state[i] = modulus(&(&state[i] + c), &self.constants.r);
        }

        new_state
    }

    pub fn cubic(&self, a: &BigInt) -> BigInt {
        modulus(&(a * a * a * a * a), &self.constants.r)
    }

    pub fn sbox(&self, state: &Vec<BigInt>, i: usize) -> Vec<BigInt> {
        let mut new_state: Vec<BigInt> = state.clone();
        if i < NROUNDSF / 2 || i >= NROUNDSF / 2 + NROUNDSP {
            for j in 0..T {
                new_state[j] = self.cubic(&state[j]);
            }
        } else {
            new_state[0] = self.cubic(&state[0]);
        }
        new_state
    }

    pub fn mix(&self, state: &Vec<BigInt>, m: &Vec<Vec<BigInt>>) -> Vec<BigInt> {
        let mut new_state: Vec<BigInt> = Vec::new();
        for i in 0..state.len() {
            new_state.push(Zero::zero());
            for j in 0..state.len() {
                new_state[i] = modulus(
                    &(&new_state[i] + modulus(&(&m[i][j] * &state[j]), &self.constants.r)),
                    &self.constants.r,
                )
            }
        }
        new_state
    }

    pub fn poseidon_hash(&self, inp: Vec<BigInt>) -> Result<BigInt, String> {
        if inp.len() == 0 || inp.len() > T {
            return Err("Wrong inputs length".to_string());
        }
        // check if arr elements are inside the finite field over R
        if !check_bigint_array_in_field(&inp, &self.constants.r) {
            return Err("elements not inside the finite field over R".to_string());
        }

        let mut state = inp.clone();
        for _ in inp.len()..T {
            state.push(Zero::zero());
        }

        for i in 0..(NROUNDSF + NROUNDSP) {
            state = self.ark(&state, &self.constants.c[i]);
            state = self.sbox(&state, i);
            state = self.mix(&state, &self.constants.m);
        }

        Ok(state[0].clone())
    }

    pub fn hash(&self, inp: Vec<BigInt>) -> Result<BigInt, String> {
        // check if arr elements are inside the finite field over R
        if !check_bigint_array_in_field(&inp, &self.constants.r) {
            return Err("elements not inside the finite field over R".to_string());
        }
        let mut r: BigInt = Zero::zero();
        for i in (0..inp.len()).step_by(5) {
            let mut five_elems: Vec<BigInt> = Vec::new();
            for j in 0..5 {
                if i + j < inp.len() {
                    five_elems.push(inp[i + j].clone());
                } else {
                    five_elems.push(Zero::zero());
                }
            }
            let ph = &self.poseidon_hash(five_elems);
            match ph {
                Result::Err(err) => return Err(err.to_string()),
                Result::Ok(res) => {
                    r = modulus(&(r + res), &self.constants.r);
                }
            }
        }
        Ok(r)
    }

    pub fn hash_bytes(&self, b: Vec<u8>) -> Result<BigInt, String> {
        let n = 31;
        let mut ints: Vec<BigInt> = Vec::new();
        for i in 0..b.len() / n {
            let v: BigInt = BigInt::from_bytes_le(Sign::Plus, &b[n * i..n * (i + 1)]);
            ints.push(v);
        }
        if b.len() % n != 0 {
            let v: BigInt = BigInt::from_bytes_le(Sign::Plus, &b[(b.len() / n) * n..]);
            ints.push(v);
        }
        self.hash(ints)
    }
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

    #[test]
    fn test_poseidon_hash() {
        let b1: BigInt = BigInt::parse_bytes(b"1", 10).unwrap();
        let b2: BigInt = BigInt::parse_bytes(b"2", 10).unwrap();
        let mut big_arr: Vec<BigInt> = Vec::new();
        big_arr.push(b1.clone());
        big_arr.push(b2.clone());
        let poseidon = Poseidon::new();
        let h = poseidon.poseidon_hash(big_arr).unwrap();
        assert_eq!(
            h.to_string(),
            "12242166908188651009877250812424843524687801523336557272219921456462821518061"
        );

        let b3: BigInt = BigInt::parse_bytes(b"3", 10).unwrap();
        let b4: BigInt = BigInt::parse_bytes(b"4", 10).unwrap();
        let mut big_arr34: Vec<BigInt> = Vec::new();
        big_arr34.push(b3.clone());
        big_arr34.push(b4.clone());
        let h34 = poseidon.poseidon_hash(big_arr34).unwrap();
        assert_eq!(
            h34.to_string(),
            "17185195740979599334254027721507328033796809509313949281114643312710535000993"
        );
    }

    #[test]
    fn test_hash_bytes() {
        let msg = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.";
        let poseidon = Poseidon::new();
        let h = poseidon.hash_bytes(msg.as_bytes().to_vec()).unwrap();
        assert_eq!(
            h.to_string(),
            "11821124228916291136371255062457365369197326845706357273715164664419275913793"
        );

        let msg2 = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum. Lorem ipsum dolor sit amet.";
        let h2 = poseidon.hash_bytes(msg2.as_bytes().to_vec()).unwrap();
        assert_eq!(
            h2.to_string(),
            "10747013384255785702102976082726575658403084163954725275481577373644732938016"
        );
    }
}
