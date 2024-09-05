use std::io::Write;

use crate::{kzg::KZG, poly_util::Polynomial};
use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use crypto_hash::{Algorithm::SHA256, Hasher};

pub struct ProdProof<E: Pairing> {
    pub com_t: E::G1,
    pub com_q: E::G1,
    pub t_wk_1: E::ScalarField,
    pub pf_t_wk_1: E::G1,
    pub t_r: E::ScalarField,
    pub pf_t_r: E::G1,
    pub t_wr: E::ScalarField,
    pub pf_t_wr: E::G1,
    pub q_r: E::ScalarField,
    pub pf_q_r: E::G1,
    pub f_wr: E::ScalarField,
    pub pf_f_wr: E::G1,
}

pub fn gen_proof<E: Pairing>(
    f: Polynomial<E::ScalarField>,
    w: E::ScalarField,
    degree_w: usize,
    com_f: E::G1,
    poly_commit: &KZG<E>,
) -> Result<ProdProof<E>, String> {
    // construct t(x), where t(1) = f(1), t(w^s) = f(w^1)*...* f(w^s) for s = 1...k-1
    let f_1 = f.eval(E::ScalarField::from(1u8));
    let mut t_evals = vec![(E::ScalarField::from(1u8), f_1)];
    // calculate t(w^i)
    for i in 1..degree_w {
        let wi = w.pow(&[i.try_into().unwrap()]);
        let t_wi = f.eval(wi) * t_evals[i - 1].1;
        t_evals.push((wi, t_wi))
    }
    let t = Polynomial::interpolate(&t_evals);
    // construct t1(X) = t(w * X) - t(X) * f(w * X)
    //t(wX)
    let mut t_wx = t.clone();
    for i in 1..t_wx.degree() + 1 {
        t_wx[i] = t_wx[i] * w.pow(&[i.try_into().unwrap()])
    }
    //f(wX)
    let mut f_wx = f.clone();
    for i in 1..f_wx.degree() + 1 {
        f_wx[i] = f_wx[i] * w.pow(&[i.try_into().unwrap()])
    }
    let t1 = &t_wx - &(&t * &f_wx);
    // vanishing polynimial of W = {1, w, w^2, ..., w^k-1} where w is kth root of unity
    let mut z_w = vec![E::ScalarField::ZERO; degree_w + 1];
    z_w[0] = -E::ScalarField::from(1u32);
    z_w[degree_w] = E::ScalarField::from(1u32);
    let z_w = Polynomial::new(z_w);

    // calculate q(x) = t1(x)/z(x)
    let (q, _) = (&t1 / &z_w)?;
    let com_t = poly_commit.commit(&t)?;
    let com_q = poly_commit.commit(&q)?;

    //public coin
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    let mut count = 0usize;
    let r;
    loop {
        match E::ScalarField::from_random_bytes(&hasher.finish()) {
            Some(coin) => {
                r = coin;
                break;
            }
            None => {
                hasher.write_all(&count.to_be_bytes()).unwrap();
                count += 1
            }
        };
    }

    // t(w^k-1)
    let (t_wk_1, pf_t_wk_1) = poly_commit.eval(&t, t_evals[degree_w - 1].0)?;
    // t(r)
    let (t_r, pf_t_r) = poly_commit.eval(&t, r)?;
    // t(wr)
    let (t_wr, pf_t_wr) = poly_commit.eval(&t, w * r)?;
    // q(r)
    let (q_r, pf_q_r) = poly_commit.eval(&q, r)?;
    // f(wr)
    let (f_wr, pf_f_wr) = poly_commit.eval(&f, w * r)?;

    Ok(ProdProof {
        com_t,
        com_q,
        t_wk_1,
        pf_t_wk_1,
        t_r,
        pf_t_r,
        t_wr,
        pf_t_wr,
        q_r,
        pf_q_r,
        f_wr,
        pf_f_wr,
    })
}

pub fn verify<E: Pairing>(
    proof: ProdProof<E>,
    w: E::ScalarField,
    degree_w: usize,
    poly_commit: &KZG<E>,
    com_f: E::G1,
) -> bool {
    //public coin
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    let mut count = 0usize;
    let r;
    loop {
        match E::ScalarField::from_random_bytes(&hasher.finish()) {
            Some(coin) => {
                r = coin;
                break;
            }
            None => {
                hasher.write_all(&count.to_be_bytes()).unwrap();
                count += 1
            }
        };
    }

    // vanishing polynimial of W = {1, w, w^2, ..., w^k-1} where w is kth root of unity
    let mut z_w = vec![E::ScalarField::ZERO; degree_w + 1];
    z_w[0] = -E::ScalarField::from(1u32);
    z_w[degree_w] = E::ScalarField::from(1u32);
    let z_w = Polynomial::new(z_w);

    let wk_1 = w.pow(&[(degree_w - 1).try_into().unwrap()]);
    // verify t(w^k-1)
    if proof.t_wk_1 != E::ScalarField::from(1u8) {
        return false;
    };
    let res = poly_commit.verify(proof.com_t, wk_1, proof.t_wk_1, proof.pf_t_wk_1);
    if !res {
        return false;
    };
    // verify t(r)
    let res = poly_commit.verify(proof.com_t, r, proof.t_r, proof.pf_t_r);
    if !res {
        return false;
    }
    // verify t(wr)
    let res = poly_commit.verify(proof.com_t, w * r, proof.t_wr, proof.pf_t_wr);
    if !res {
        return false;
    }
    // verify q(r)
    let res = poly_commit.verify(proof.com_q, r, proof.q_r, proof.pf_q_r);
    if !res {
        return false;
    }
    // verify f(wr)
    let res = poly_commit.verify(com_f, w * r, proof.f_wr, proof.pf_f_wr);
    if !res {
        return false;
    }
    // check if t(wr) - t(r)f(wr) = q(r)*(r^k-1)
    if proof.t_wr - (proof.t_r * proof.f_wr) != proof.q_r * z_w.eval(r) {
        return false;
    }

    true
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

        let f = Polynomial::interpolate(&[
            (w.pow(&[1]), Fr::from(2)),
            (w.pow(&[2]), Fr::from(3)),
            (w.pow(&[3]), Fr::from(4)),
            (w.pow(&[4]), Fr::from(5)),
            (w.pow(&[5]), Fr::from(120).inverse().unwrap()),
            (w.pow(&[6]), Fr::from(7)),
            (w.pow(&[7]), Fr::from(8)),
            (w.pow(&[8]), Fr::from(56).inverse().unwrap()),
            (Fr::from(10), Fr::from(12)),
        ]);
        let mut poly_commit =
            KZG::<Bls12_381>::new(G1Projective::generator(), G2Projective::generator(), 8);
        poly_commit.random_setup(&mut rand::thread_rng());
        let com_f = poly_commit.commit(&f).unwrap();

        let pf = gen_proof(f, w, degree_w, com_f, &poly_commit).unwrap();
        let res = verify(pf, w, degree_w, &poly_commit, com_f);
        assert!(res);
    }

    #[test]
    pub fn test_prod_soundness() {
        // create w as an 8th root of unity
        let modulus: BigUint = Fr::MODULUS.into();
        let one: num_bigint::BigUint = Fr::one().into();
        let ord = &modulus - one;
        let fr_gen = Fr::from(7);
        let subgroup_ord = Fr::from(8);
        let w = fr_gen.pow((ord / BigUint::from(subgroup_ord)).to_u64_digits());
        let degree_w = 8;

        let f = Polynomial::interpolate(&[
            (w.pow(&[1]), Fr::from(10)),
            (w.pow(&[2]), Fr::from(3)),
            (w.pow(&[3]), Fr::from(4)),
            (w.pow(&[4]), Fr::from(5)),
            (w.pow(&[5]), Fr::from(120).inverse().unwrap()),
            (w.pow(&[6]), Fr::from(7)),
            (w.pow(&[7]), Fr::from(8)),
            (w.pow(&[8]), Fr::from(56).inverse().unwrap()),
            (Fr::from(10), Fr::from(12)),
        ]);
        let mut poly_commit =
            KZG::<Bls12_381>::new(G1Projective::generator(), G2Projective::generator(), 8);
        poly_commit.random_setup(&mut rand::thread_rng());
        let com_f = poly_commit.commit(&f).unwrap();

        let pf = gen_proof(f, w, degree_w, com_f, &poly_commit).unwrap();
        let res = verify(pf, w, degree_w, &poly_commit, com_f);
        assert!(!res);
    }
}
