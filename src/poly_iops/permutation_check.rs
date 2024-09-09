use crate::{kzg::KZG, poly_util::Polynomial};
use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use crypto_hash::{Algorithm::SHA256, Hasher};
use std::io::Write;

use super::prod_check;

pub struct PerProof<E: Pairing> {
    pub com_f1: E::G1,
    pub com_g1: E::G1,
    pub f1_rf1: E::ScalarField,
    pub pf_f1_rf1: E::G1,
    pub per_rf1: E::ScalarField,
    pub pf_per_rf1: E::G1,
    pub f_rf1: E::ScalarField,
    pub pf_f_rf1: E::G1,
    pub g1_rg1: E::ScalarField,
    pub pf_g1_rg1: E::G1,
    pub g_rg1: E::ScalarField,
    pub pf_g_rg1: E::G1,
    pub prod_pf: prod_check::RationalFuncProdProof<E>,
}

// proof that g(W) is the same as f(W), permuted by the prescribed per function
// where per is a permutation on W = {1, w^1, ... w^degree_w-1}
pub fn gen_proof<E: Pairing>(
    f: Polynomial<E::ScalarField>,
    g: Polynomial<E::ScalarField>,
    per: Polynomial<E::ScalarField>,
    w: E::ScalarField,
    degree_w: usize,
    com_f: E::G1,
    com_g: E::G1,
    com_per: E::G1,
    poly_commit: &KZG<E>,
) -> Result<PerProof<E>, String> {
    // generate public coin r_per, s_per, r_f1, r_g1
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    com_g.serialize_uncompressed(&mut hasher).unwrap();
    com_per.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    hasher.write_all(b"R_PER").unwrap();
    let mut count = 0usize;
    let r_per;
    loop {
        match E::ScalarField::from_random_bytes(&hasher.finish()) {
            Some(coin) => {
                r_per = coin;
                break;
            }
            None => {
                hasher.write_all(&count.to_be_bytes()).unwrap();
                count += 1
            }
        };
    }
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    com_g.serialize_uncompressed(&mut hasher).unwrap();
    com_per.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    hasher.write_all(b"S_PER").unwrap();
    let mut count = 0usize;
    let s_per;
    loop {
        match E::ScalarField::from_random_bytes(&hasher.finish()) {
            Some(coin) => {
                s_per = coin;
                break;
            }
            None => {
                hasher.write_all(&count.to_be_bytes()).unwrap();
                count += 1
            }
        };
    }
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    com_g.serialize_uncompressed(&mut hasher).unwrap();
    com_per.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    hasher.write_all(b"R_F1").unwrap();
    let mut count = 0usize;
    let r_f1;
    loop {
        match E::ScalarField::from_random_bytes(&hasher.finish()) {
            Some(coin) => {
                r_f1 = coin;
                break;
            }
            None => {
                hasher.write_all(&count.to_be_bytes()).unwrap();
                count += 1
            }
        };
    }
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    com_g.serialize_uncompressed(&mut hasher).unwrap();
    com_per.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    hasher.write_all(b"R_G1").unwrap();
    let mut count = 0usize;
    let r_g1;
    loop {
        match E::ScalarField::from_random_bytes(&hasher.finish()) {
            Some(coin) => {
                r_g1 = coin;
                break;
            }
            None => {
                hasher.write_all(&count.to_be_bytes()).unwrap();
                count += 1
            }
        };
    }

    // f1 = (s * per(x) + f(x) - r) & g1 = (sx + g(x) - r)
    let s = Polynomial::new(vec![s_per]);
    let sx = Polynomial::new(vec![E::ScalarField::ZERO, s_per]);
    let r = Polynomial::new(vec![r_per]);
    let f1 = &(&(&s * &per) + &f) - &r;
    let g1 = &(&sx + &g) - &r;
    let com_f1 = poly_commit.commit(&f1)?;
    let com_g1 = poly_commit.commit(&g1)?;
    // prove that f1 & g1 is constructed correctly
    // f1(r_f1)
    let (f1_rf1, pf_f1_rf1) = poly_commit.eval(&f1, r_f1)?;
    // per(r_f1)
    let (per_rf1, pf_per_rf1) = poly_commit.eval(&per, r_f1)?;
    // f(r_f1)
    let (f_rf1, pf_f_rf1) = poly_commit.eval(&f, r_f1)?;
    // g1(r_g1)
    let (g1_rg1, pf_g1_rg1) = poly_commit.eval(&g1, r_g1)?;
    // g(r_g1)
    let (g_rg1, pf_g_rg1) = poly_commit.eval(&g, r_g1)?;

    // prod proof f1 / g1
    let prod_pf =
        prod_check::rational_func_gen_proof(f1, g1, w, degree_w, com_f1, com_g1, poly_commit)?;

    Ok(PerProof {
        com_f1,
        com_g1,
        f1_rf1,
        pf_f1_rf1,
        per_rf1,
        pf_per_rf1,
        f_rf1,
        pf_f_rf1,
        g1_rg1,
        pf_g1_rg1,
        g_rg1,
        pf_g_rg1,
        prod_pf,
    })
}

pub fn verify<E: Pairing>(
    pf: PerProof<E>,
    w: E::ScalarField,
    degree_w: usize,
    com_f: E::G1,
    com_g: E::G1,
    com_per: E::G1,
    poly_commit: &KZG<E>,
) -> bool {
    // generate public coin r_per, s_per, r_f1, r_g1
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    com_g.serialize_uncompressed(&mut hasher).unwrap();
    com_per.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    hasher.write_all(b"R_PER").unwrap();
    let mut count = 0usize;
    let r_per;
    loop {
        match E::ScalarField::from_random_bytes(&hasher.finish()) {
            Some(coin) => {
                r_per = coin;
                break;
            }
            None => {
                hasher.write_all(&count.to_be_bytes()).unwrap();
                count += 1
            }
        };
    }
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    com_g.serialize_uncompressed(&mut hasher).unwrap();
    com_per.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    hasher.write_all(b"S_PER").unwrap();
    let mut count = 0usize;
    let s_per;
    loop {
        match E::ScalarField::from_random_bytes(&hasher.finish()) {
            Some(coin) => {
                s_per = coin;
                break;
            }
            None => {
                hasher.write_all(&count.to_be_bytes()).unwrap();
                count += 1
            }
        };
    }
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    com_g.serialize_uncompressed(&mut hasher).unwrap();
    com_per.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    hasher.write_all(b"R_F1").unwrap();
    let mut count = 0usize;
    let r_f1;
    loop {
        match E::ScalarField::from_random_bytes(&hasher.finish()) {
            Some(coin) => {
                r_f1 = coin;
                break;
            }
            None => {
                hasher.write_all(&count.to_be_bytes()).unwrap();
                count += 1
            }
        };
    }
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    com_g.serialize_uncompressed(&mut hasher).unwrap();
    com_per.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    hasher.write_all(b"R_G1").unwrap();
    let mut count = 0usize;
    let r_g1;
    loop {
        match E::ScalarField::from_random_bytes(&hasher.finish()) {
            Some(coin) => {
                r_g1 = coin;
                break;
            }
            None => {
                hasher.write_all(&count.to_be_bytes()).unwrap();
                count += 1
            }
        };
    }

    let res = prod_check::rational_func_verify(
        pf.prod_pf,
        w,
        degree_w,
        poly_commit,
        pf.com_f1,
        pf.com_g1,
    );
    if !res {
        return false;
    }
    let res = poly_commit.verify(pf.com_f1, r_f1, pf.f1_rf1, pf.pf_f1_rf1);
    if !res {
        return false;
    }
    let res = poly_commit.verify(com_per, r_f1, pf.per_rf1, pf.pf_per_rf1);
    if !res {
        return false;
    }
    let res = poly_commit.verify(com_f, r_f1, pf.f_rf1, pf.pf_f_rf1);
    if !res {
        return false;
    }
    let res = poly_commit.verify(pf.com_g1, r_g1, pf.g1_rg1, pf.pf_g1_rg1);
    if !res {
        return false;
    }
    let res = poly_commit.verify(com_g, r_g1, pf.g_rg1, pf.pf_g_rg1);
    if !res {
        return false;
    }
    // check f1(x) = (s * per(x) + f(x) - r) & g1(x) = (sx + g(x) - r)
    if !(pf.f1_rf1 == s_per * pf.per_rf1 + pf.f_rf1 - r_per) {
        return false;
    }
    if !(pf.g1_rg1 == (s_per * r_g1) + pf.g_rg1 - r_per) {
        return false;
    }

    return true;
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bls12_381::{Bls12_381, Fr, G1Projective, G2Projective};
    use ark_ec::Group;
    use ark_ff::PrimeField;
    use ark_std::One;
    use num_bigint::BigUint;

    #[test]
    pub fn happy_path() {
        // create w as an 8th root of unity
        let modulus: BigUint = Fr::MODULUS.into();
        let one: num_bigint::BigUint = Fr::one().into();
        let ord = &modulus - one;
        let fr_gen = Fr::from(7);
        let subgroup_ord = Fr::from(8);
        let w = fr_gen.pow((ord / BigUint::from(subgroup_ord)).to_u64_digits());
        let degree_w = 8;

        // create permutation function
        // per(w^0) = w^2, per(w^1) = w^5, per(w^2) = w^6, per(w^3) = w^0
        // per(w^4) = w^1, per(w^5) = w^3, per(w^6) = w^7, per(w^7) = w^4
        let per = Polynomial::interpolate(&[
            (w.pow(&[0]), w.pow(&[2])),
            (w.pow(&[1]), w.pow(&[5])),
            (w.pow(&[2]), w.pow(&[6])),
            (w.pow(&[3]), w.pow(&[0])),
            (w.pow(&[4]), w.pow(&[1])),
            (w.pow(&[5]), w.pow(&[3])),
            (w.pow(&[6]), w.pow(&[7])),
            (w.pow(&[7]), w.pow(&[4])),
        ]);

        let f = Polynomial::interpolate(&[
            (w.pow(&[0]), Fr::from(12)),
            (w.pow(&[1]), Fr::from(40)),
            (w.pow(&[2]), Fr::from(2)),
            (w.pow(&[3]), Fr::from(9)),
            (w.pow(&[4]), Fr::from(300).inverse().unwrap()),
            (w.pow(&[5]), Fr::from(24)),
            (w.pow(&[6]), Fr::from(1)),
            (w.pow(&[7]), Fr::from(17).inverse().unwrap()),
            (Fr::from(10), Fr::from(25)),
        ]);

        let g = Polynomial::interpolate(&[
            (w.pow(&[2]), Fr::from(12)),
            (w.pow(&[5]), Fr::from(40)),
            (w.pow(&[6]), Fr::from(2)),
            (w.pow(&[0]), Fr::from(9)),
            (w.pow(&[1]), Fr::from(300).inverse().unwrap()),
            (w.pow(&[3]), Fr::from(24)),
            (w.pow(&[7]), Fr::from(1)),
            (w.pow(&[4]), Fr::from(17).inverse().unwrap()),
            (Fr::from(17), Fr::from(35)),
        ]);

        let mut poly_commit =
            KZG::<Bls12_381>::new(G1Projective::generator(), G2Projective::generator(), 8);
        poly_commit.random_setup(&mut rand::thread_rng());
        let com_f = poly_commit.commit(&f).unwrap();
        let com_g = poly_commit.commit(&g).unwrap();
        let com_per = poly_commit.commit(&per).unwrap();

        let pf = gen_proof(f, g, per, w, degree_w, com_f, com_g, com_per, &poly_commit).unwrap();
        let res = verify(pf, w, degree_w, com_f, com_g, com_per, &poly_commit);
        assert!(res);
    }

    #[test]
    pub fn incorrect_per_func() {
        // create w as an 8th root of unity
        let modulus: BigUint = Fr::MODULUS.into();
        let one: num_bigint::BigUint = Fr::one().into();
        let ord = &modulus - one;
        let fr_gen = Fr::from(7);
        let subgroup_ord = Fr::from(8);
        let w = fr_gen.pow((ord / BigUint::from(subgroup_ord)).to_u64_digits());
        let degree_w = 8;

        // create permutation function
        // wrong permutation
        // per(w^0) = w^2, per(w^1) = w^6, per(w^2) = w^5, per(w^3) = w^0
        // per(w^4) = w^1, per(w^5) = w^3, per(w^6) = w^7, per(w^7) = w^4
        let per = Polynomial::interpolate(&[
            (w.pow(&[0]), w.pow(&[2])),
            (w.pow(&[1]), w.pow(&[6])),
            (w.pow(&[2]), w.pow(&[5])),
            (w.pow(&[3]), w.pow(&[0])),
            (w.pow(&[4]), w.pow(&[1])),
            (w.pow(&[5]), w.pow(&[3])),
            (w.pow(&[6]), w.pow(&[7])),
            (w.pow(&[7]), w.pow(&[4])),
        ]);

        let f = Polynomial::interpolate(&[
            (w.pow(&[0]), Fr::from(12)),
            (w.pow(&[1]), Fr::from(40)),
            (w.pow(&[2]), Fr::from(2)),
            (w.pow(&[3]), Fr::from(9)),
            (w.pow(&[4]), Fr::from(300).inverse().unwrap()),
            (w.pow(&[5]), Fr::from(24)),
            (w.pow(&[6]), Fr::from(1)),
            (w.pow(&[7]), Fr::from(17).inverse().unwrap()),
            (Fr::from(10), Fr::from(25)),
        ]);

        let g = Polynomial::interpolate(&[
            (w.pow(&[2]), Fr::from(12)),
            (w.pow(&[5]), Fr::from(40)),
            (w.pow(&[6]), Fr::from(2)),
            (w.pow(&[0]), Fr::from(9)),
            (w.pow(&[1]), Fr::from(300).inverse().unwrap()),
            (w.pow(&[3]), Fr::from(24)),
            (w.pow(&[7]), Fr::from(1)),
            (w.pow(&[4]), Fr::from(17).inverse().unwrap()),
            (Fr::from(17), Fr::from(35)),
        ]);

        let mut poly_commit =
            KZG::<Bls12_381>::new(G1Projective::generator(), G2Projective::generator(), 8);
        poly_commit.random_setup(&mut rand::thread_rng());
        let com_f = poly_commit.commit(&f).unwrap();
        let com_g = poly_commit.commit(&g).unwrap();
        let com_per = poly_commit.commit(&per).unwrap();

        let pf = gen_proof(f, g, per, w, degree_w, com_f, com_g, com_per, &poly_commit).unwrap();
        let res = verify(pf, w, degree_w, com_f, com_g, com_per, &poly_commit);
        assert!(!res);
    }

    #[test]
    pub fn incorrect_g_func() {
        // create w as an 8th root of unity
        let modulus: BigUint = Fr::MODULUS.into();
        let one: num_bigint::BigUint = Fr::one().into();
        let ord = &modulus - one;
        let fr_gen = Fr::from(7);
        let subgroup_ord = Fr::from(8);
        let w = fr_gen.pow((ord / BigUint::from(subgroup_ord)).to_u64_digits());
        let degree_w = 8;

        // create permutation function
        // wrong permutation
        // per(w^0) = w^2, per(w^1) = w^6, per(w^2) = w^5, per(w^3) = w^0
        // per(w^4) = w^1, per(w^5) = w^3, per(w^6) = w^7, per(w^7) = w^4
        let per = Polynomial::interpolate(&[
            (w.pow(&[0]), w.pow(&[2])),
            (w.pow(&[1]), w.pow(&[6])),
            (w.pow(&[2]), w.pow(&[5])),
            (w.pow(&[3]), w.pow(&[0])),
            (w.pow(&[4]), w.pow(&[1])),
            (w.pow(&[5]), w.pow(&[3])),
            (w.pow(&[6]), w.pow(&[7])),
            (w.pow(&[7]), w.pow(&[4])),
        ]);

        let f = Polynomial::interpolate(&[
            (w.pow(&[0]), Fr::from(12)),
            (w.pow(&[1]), Fr::from(40)),
            (w.pow(&[2]), Fr::from(2)),
            (w.pow(&[3]), Fr::from(9)),
            (w.pow(&[4]), Fr::from(300).inverse().unwrap()),
            (w.pow(&[5]), Fr::from(24)),
            (w.pow(&[6]), Fr::from(1)),
            (w.pow(&[7]), Fr::from(17).inverse().unwrap()),
            (Fr::from(10), Fr::from(25)),
        ]);

        let g = Polynomial::interpolate(&[
            (w.pow(&[2]), Fr::from(12)),
            (w.pow(&[5]), Fr::from(40)),
            (w.pow(&[6]), Fr::from(1)),
            (w.pow(&[0]), Fr::from(9)),
            (w.pow(&[1]), Fr::from(300).inverse().unwrap()),
            (w.pow(&[3]), Fr::from(24)),
            (w.pow(&[7]), Fr::from(1)),
            (w.pow(&[4]), Fr::from(17).inverse().unwrap()),
            (Fr::from(17), Fr::from(35)),
        ]);

        let mut poly_commit =
            KZG::<Bls12_381>::new(G1Projective::generator(), G2Projective::generator(), 8);
        poly_commit.random_setup(&mut rand::thread_rng());
        let com_f = poly_commit.commit(&f).unwrap();
        let com_g = poly_commit.commit(&g).unwrap();
        let com_per = poly_commit.commit(&per).unwrap();

        let pf = gen_proof(f, g, per, w, degree_w, com_f, com_g, com_per, &poly_commit).unwrap();
        let res = verify(pf, w, degree_w, com_f, com_g, com_per, &poly_commit);
        assert!(!res);
    }
}
