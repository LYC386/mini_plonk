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
