use std::io::Write;

use crate::{kzg::KZG, poly_util::Polynomial};
use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use crypto_hash::{Algorithm::SHA256, Hasher};

pub fn gen_proof<E: Pairing>(
    f: Polynomial<E::ScalarField>,
    w: E::ScalarField,
    degree_w: usize,
    com_f: E::G1,
    poly_commit: &KZG<E>,
) -> Result<
    (
        E::G1,
        E::G1,
        E::ScalarField,
        E::G1,
        E::ScalarField,
        E::G1,
        E::ScalarField,
        E::G1,
        E::ScalarField,
        E::G1,
        E::ScalarField,
        E::G1,
    ),
    String,
> {
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
    for i in 1..t_wx.degree() {
        t_wx[i] = t_wx[i] * w.pow(&[i.try_into().unwrap()])
    }
    //f(wX)
    let mut f_wx = f.clone();
    for i in 1..f_wx.degree() {
        f_wx[i] = f_wx[i] * w.pow(&[i.try_into().unwrap()])
    }
    let t1 = &t_wx - &(&t * &f_wx);

    // vanishing polynimial of W = P{1, w, w^2, ..., w^k-1}
    // TODO: set w as kth power of root
    let mut z_w = Polynomial::new(vec![E::ScalarField::ONE]);
    for i in 0..degree_w {
        z_w = &z_w * &Polynomial::new(vec![-w.pow(&[i.try_into().unwrap()]), E::ScalarField::ONE]);
    }
    // calculate q(x) = t1(x)/z(x)
    let (q, _) = (&t1 / &z_w)?;
    let com_t = poly_commit.commit(&t)?;
    let com_q = poly_commit.commit(&q)?;

    //public coin
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    let r = E::ScalarField::from_random_bytes(&hasher.finish()).unwrap();
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

    Ok((
        com_t, com_q, t_wk_1, pf_t_wk_1, t_r, pf_t_r, t_wr, pf_t_wr, q_r, pf_q_r, f_wr, pf_f_wr,
    ))
}
