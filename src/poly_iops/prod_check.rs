use ark_ec::pairing::Pairing;
use ark_ff::Field;

use crate::{kzg::KZG, poly_util::Polynomial};

// pub fn gen_proof<E: Pairing>(
//     f: Polynomial<E::ScalarField>,
//     w: E::ScalarField,
//     degree_w: usize,
//     com_f: E::G1,
//     poly_commit: &KZG<E>,
// ) -> Result<(), String> {
//     // construct t(x), where t(1) = f(1), t(w^s) = f(w^1)*...* f(w^s) for s = 1...k-1
//     let f_1 = f.eval(E::ScalarField::from(1u8));
//     let mut t_evals = vec![(E::ScalarField::from(1u8), f_1)];
//     // calculate t(w^i)
//     for i in 1..degree_w {
//         let wi = w.pow(&[i.try_into().unwrap()]);
//         let t_wi = f.eval(wi) * t_evals[i - 1].1;
//         t_evals.push((wi, t_wi))
//     }
//     let t = Polynomial::interpolate(&t_evals);
//     // construct t1(X) = t(w * X) - t(X) * f(w * X)
//     //t(wX)
//     let mut t_wx = t.clone();
//     for i in 1..t_wx.degree() {
//         t_wx[i] = t_wx[i] * w.pow(&[i.try_into().unwrap()])
//     }
//     //f(wX)
//     let mut f_wx = f.clone();
//     for i in 1..f_wx.degree() {
//         f_wx[i] = f_wx[i] * w.pow(&[i.try_into().unwrap()])
//     }
//     let t1 = &t_wx - &(&t * &f_wx);

//     // vanishing polynimial of W = P{1, w, w^2, ..., w^k-1}
//     // TODO: set w as kth power of root
//     let mut z_w = Polynomial::new(vec![E::ScalarField::ONE]);
//     for i in 0..degree_w {
//         z_w = &z_w * &Polynomial::new(vec![-w.pow(&[i.try_into().unwrap()]), E::ScalarField::ONE]);
//     }
//     // calculate q(x) = t1(x)/z(x)
//     let q = &t1 / &z_w;
//     // let com_t = poly_commit.commit(t)?;
// }
