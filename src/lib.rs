extern crate rand;
#[macro_use]
extern crate ff;
use ff::*;

#[derive(PrimeField)]
#[PrimeFieldModulus = "21888242871839275222246405745257275088548364400416034343698204186575808495617"]
#[PrimeFieldGenerator = "7"]
pub struct Fr(FrRepr);

mod constants;

#[derive(Debug)]
pub struct Constants {
    pub c: Vec<Vec<Fr>>,
    pub m: Vec<Vec<Vec<Fr>>>,
    pub n_rounds_f: usize,
    pub n_rounds_p: Vec<usize>,
}
pub fn load_constants() -> Constants {
    let (c_str, m_str) = constants::constants();
    let mut c: Vec<Vec<Fr>> = Vec::new();
    for i in 0..c_str.len() {
        let mut cci: Vec<Fr> = Vec::new();
        for j in 0..c_str[i].len() {
            let b: Fr = Fr::from_str(c_str[i][j]).unwrap();
            cci.push(b);
        }
        c.push(cci);
    }
    let mut m: Vec<Vec<Vec<Fr>>> = Vec::new();
    for i in 0..m_str.len() {
        let mut mi: Vec<Vec<Fr>> = Vec::new();
        for j in 0..m_str[i].len() {
            let mut mij: Vec<Fr> = Vec::new();
            for k in 0..m_str[i][j].len() {
                let b: Fr = Fr::from_str(m_str[i][j][k]).unwrap();
                mij.push(b);
            }
            mi.push(mij);
        }
        m.push(mi);
    }
    Constants {
        c: c,
        m: m,
        n_rounds_f: 8,
        n_rounds_p: vec![56, 57, 56, 60, 60, 63, 64, 63],
    }
}

pub struct Poseidon {
    constants: Constants,
}
impl Poseidon {
    pub fn new() -> Poseidon {
        Poseidon {
            constants: load_constants(),
        }
    }
    pub fn ark(&self, state: &mut Vec<Fr>, c: &Vec<Fr>, it: usize) {
        for i in 0..state.len() {
            state[i].add_assign(&c[it + i]);
        }
    }

    pub fn sbox(&self, n_rounds_f: usize, n_rounds_p: usize, state: &mut Vec<Fr>, i: usize) {
        if i < n_rounds_f / 2 || i >= n_rounds_f / 2 + n_rounds_p {
            for j in 0..state.len() {
                let aux = state[j];
                state[j].square();
                state[j].square();
                state[j].mul_assign(&aux);
            }
        } else {
            let aux = state[0];
            state[0].square();
            state[0].square();
            state[0].mul_assign(&aux);
        }
    }

    pub fn mix(&self, state: &Vec<Fr>, m: &Vec<Vec<Fr>>) -> Vec<Fr> {
        let mut new_state: Vec<Fr> = Vec::new();
        for i in 0..state.len() {
            new_state.push(Fr::zero());
            for j in 0..state.len() {
                let mut mij = m[j][i];
                mij.mul_assign(&state[j]);
                new_state[i].add_assign(&mij);
            }
        }
        new_state.clone()
    }

    pub fn hash(&self, inp: Vec<Fr>) -> Result<Fr, String> {
        let t = inp.len() + 1;
        if inp.len() == 0 || inp.len() >= self.constants.n_rounds_p.len() - 1 {
            return Err("Wrong inputs length".to_string());
        }
        let n_rounds_f = self.constants.n_rounds_f.clone();
        let n_rounds_p = self.constants.n_rounds_p[t - 2].clone();

        let mut state = inp.clone();
        for _ in inp.len()..t {
            state.push(Fr::zero());
        }
        // state[state.len() - 1] = Fr::zero();

        for i in 0..(n_rounds_f + n_rounds_p) {
            self.ark(&mut state, &self.constants.c[t - 2], i * t);
            self.sbox(n_rounds_f, n_rounds_p, &mut state, i);
            if i < n_rounds_f + n_rounds_p - 1 {
                state = self.mix(&state, &self.constants.m[t - 2]);
            }
        }

        Ok(state[0])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ff() {
        let a = Fr::from_repr(FrRepr::from(2)).unwrap();
        assert_eq!(
            "0000000000000000000000000000000000000000000000000000000000000002",
            to_hex(&a)
        );

        let b: Fr = Fr::from_str(
            "21888242871839275222246405745257275088548364400416034343698204186575808495619",
        )
        .unwrap();
        assert_eq!(
            "0000000000000000000000000000000000000000000000000000000000000002",
            to_hex(&b)
        );
        assert_eq!(&a, &b);
    }

    #[test]
    fn test_load_constants() {
        let cons = load_constants();
        assert_eq!(
            cons.c[0][0].to_string(),
            "Fr(0x09c46e9ec68e9bd4fe1faaba294cba38a71aa177534cdd1b6c7dc0dbd0abd7a7)"
        );
        assert_eq!(
            cons.c[cons.c.len() - 1][0].to_string(),
            "Fr(0x2088ce9534577bf38be7bc457f2756d558d66e0c07b9cc001a580bd42cda0e77)"
        );
        assert_eq!(
            cons.m[0][0][0].to_string(),
            "Fr(0x066f6f85d6f68a85ec10345351a23a3aaf07f38af8c952a7bceca70bd2af7ad5)"
        );
        assert_eq!(
            cons.m[cons.m.len() - 1][0][0].to_string(),
            "Fr(0x0190f922d97c8a7dcf0a142a3be27749d1c64bc22f1c556aaa24925d158cac56)"
        );
    }

    #[test]
    fn test_hash() {
        let b0: Fr = Fr::from_str("0").unwrap();
        let b1: Fr = Fr::from_str("1").unwrap();
        let b2: Fr = Fr::from_str("2").unwrap();
        let b3: Fr = Fr::from_str("3").unwrap();
        let b4: Fr = Fr::from_str("4").unwrap();
        let b5: Fr = Fr::from_str("5").unwrap();
        let b6: Fr = Fr::from_str("6").unwrap();

        let mut big_arr: Vec<Fr> = Vec::new();
        big_arr.push(b1.clone());
        let poseidon = Poseidon::new();
        let h = poseidon.hash(big_arr.clone()).unwrap();
        assert_eq!(
            h.to_string(),
            "Fr(0x186a5454a7c47c73dfc74ac32ea40a57d27eeb4e2bfc6551dd7b66686d3fd1ab)" // "11043376183861534927536506085090418075369306574649619885724436265926427398571"
        );

        let mut big_arr: Vec<Fr> = Vec::new();
        big_arr.push(b1.clone());
        big_arr.push(b2.clone());
        let poseidon = Poseidon::new();
        let h = poseidon.hash(big_arr.clone()).unwrap();
        assert_eq!(
            h.to_string(),
            "Fr(0x25d86fb7c42fd70a7e800e871f22f2f03a282abb18f86c347a1078a92f713f60)" // "17117985411748610629288516079940078114952304104811071254131751175361957805920"
        );

        let mut big_arr: Vec<Fr> = Vec::new();
        big_arr.push(b1.clone());
        big_arr.push(b2.clone());
        big_arr.push(b0.clone());
        big_arr.push(b0.clone());
        big_arr.push(b0.clone());
        let poseidon = Poseidon::new();
        let h = poseidon.hash(big_arr.clone()).unwrap();
        assert_eq!(
            h.to_string(),
            "Fr(0x08ca0a9154fccd6426092b2404e1ceeb80a7849734f1d3fe7952c2075e489566)" // "3975478831357328722254985704342968745327876719981393787143845259590563829094"
        );

        let mut big_arr: Vec<Fr> = Vec::new();
        big_arr.push(b1.clone());
        big_arr.push(b2.clone());
        big_arr.push(b0.clone());
        big_arr.push(b0.clone());
        big_arr.push(b0.clone());
        big_arr.push(b0.clone());
        let poseidon = Poseidon::new();
        let h = poseidon.hash(big_arr.clone()).unwrap();
        assert_eq!(
            h.to_string(),
            "Fr(0x2bb6c270db4ca49d129e315cdad9e0e678c1692c420dbf4667fdabc0f158e4ae)" // "19772360636270345724087386688434825760738403416279047262510528378903625000110"
        );

        let mut big_arr: Vec<Fr> = Vec::new();
        big_arr.push(b3.clone());
        big_arr.push(b4.clone());
        big_arr.push(b0.clone());
        big_arr.push(b0.clone());
        big_arr.push(b0.clone());
        let poseidon = Poseidon::new();
        let h = poseidon.hash(big_arr.clone()).unwrap();
        assert_eq!(
            h.to_string(),
            "Fr(0x07087ef123b0fc18a7487a9b3112aec23601e3d2b7ea27a85b35c7ecb595e6f6)" // "3181200837746671699652342497997860344148947482942465819251904554707352676086"
        );

        let mut big_arr: Vec<Fr> = Vec::new();
        big_arr.push(b3.clone());
        big_arr.push(b4.clone());
        big_arr.push(b0.clone());
        big_arr.push(b0.clone());
        big_arr.push(b0.clone());
        big_arr.push(b0.clone());
        let h = poseidon.hash(big_arr.clone()).unwrap();
        assert_eq!(
            h.to_string(),
            "Fr(0x128a815839bb66db834533b9c837e5a09df55e90aa9aba7ad46782234e083c20)" // "8386348873272147968934270337233829407378789978142456170950021426339096575008"
        );

        let mut big_arr: Vec<Fr> = Vec::new();
        big_arr.push(b1.clone());
        big_arr.push(b2.clone());
        big_arr.push(b3.clone());
        big_arr.push(b4.clone());
        big_arr.push(b5.clone());
        big_arr.push(b6.clone());
        let h = poseidon.hash(big_arr.clone()).unwrap();
        assert_eq!(
            h.to_string(),
            "Fr(0x0b807dafd5ecc62acdf7ae48e3a1dfb14ccc1ce398f865ac85ff0b4afd90ea6c)" // "5202465217520500374834597824465244016759843635092906214933648999760272616044"
        );
    }
}
