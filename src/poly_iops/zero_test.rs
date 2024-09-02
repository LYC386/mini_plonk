use crate::{kzg::KZG, poly_util::Polynomial};
use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use crypto_hash::{Algorithm::SHA256, Hasher};
use std::io::Write;

// generate a zero proof, returns com_q, q(r), pi_q, f(r), pi_f
pub fn gen_proof<E: Pairing>(
    w: E::ScalarField,
    f: &Polynomial<E::ScalarField>,
    com_f: E::G1,
    degree_w: usize,
    poly_commit: &KZG<E>,
) -> Result<(E::G1, E::ScalarField, E::G1, E::ScalarField, E::G1), String> {
    // vanishing polynimial of W = P{1, w, w^2, ..., w^k-1}
    // TODO: set w as kth power of root
    // let mut z_w: Vec<<E as Pairing>::ScalarField> = vec![E::ScalarField::ZERO; degree_w];
    // z_w[0] = E::ScalarField::from(1u32).neg();
    // z_w[degree_w - 1] = E::ScalarField::from(1u32);
    let mut z_w = Polynomial::new(vec![E::ScalarField::ONE]);
    for i in 0..degree_w {
        z_w = &z_w * &Polynomial::new(vec![-w.pow(&[i.try_into().unwrap()]), E::ScalarField::ONE]);
    }

    // calaulate q(x)
    let (q, _) = (f / &z_w)?;
    let com_q = poly_commit.commit(&q)?;

    // generate r
    //TODO - check what should be in the public coin
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    // for i in f {
    //     i.serialize_uncompressed(hasher).unwrap();
    // }
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    let r = E::ScalarField::from_random_bytes(&hasher.finish()).unwrap();

    // eval q(r), f(r)
    let (fr, pi_fr) = poly_commit.eval(f, r)?;
    let (qr, pi_qr) = poly_commit.eval(&q, r)?;

    return Ok((com_q, qr, pi_qr, fr, pi_fr));
}

// verify a zero proof
pub fn verify<E: Pairing>(
    w: E::ScalarField,
    com_f: E::G1,
    com_q: E::G1,
    qr: E::ScalarField,
    pi_qr: E::G1,
    fr: E::ScalarField,
    pi_fr: E::G1,
    degree_w: usize,
    poly_commit: &KZG<E>,
) -> bool {
    // vanishing polynimial of W = P{1, w, w^2, ..., w^k-1}
    // TODO: set w as kth power of root
    // let mut z_w: Vec<<E as Pairing>::ScalarField> = vec![E::ScalarField::ZERO; degree_w];
    // z_w[0] = E::ScalarField::from(1u32).neg();
    // z_w[degree_w - 1] = E::ScalarField::from(1u32);
    let mut z_w = Polynomial::new(vec![E::ScalarField::ONE]);
    for i in 0..degree_w {
        z_w = &z_w * &Polynomial::new(vec![-w.pow(&[i.try_into().unwrap()]), E::ScalarField::ONE]);
    }

    // calculate coin r
    //TODO - check what should be in the public coin
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    // for i in f {
    //     i.serialize_uncompressed(hasher).unwrap();
    // }
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    let r = E::ScalarField::from_random_bytes(&hasher.finish()).unwrap();

    //verify qr, fr
    if !poly_commit.verify(com_f, r, fr, pi_fr) {
        return false;
    }
    if !poly_commit.verify(com_q, r, qr, pi_qr) {
        return false;
    }

    //check if f(r)=q(r)*z_w(r)
    let z_wr = z_w.eval(r);
    if fr != qr * z_wr {
        return false;
    }
    true
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bls12_381::{Bls12_381, Fr, G1Projective, G2Projective};
    use ark_ec::Group;

    #[test]
    pub fn zero_test_happy_path() {
        // W = (1,5,25)
        let w = Fr::from(5u8);
        let degree_w = 3;
        // f = (x-1)(x-5)(x-25)(x+3)
        let f = Polynomial::new(vec![
            Fr::from(-375),
            Fr::from(340),
            Fr::from(62),
            Fr::from(-28),
            Fr::from(1),
        ]);
        let mut poly_commit =
            KZG::<Bls12_381>::new(G1Projective::generator(), G2Projective::generator(), 4);
        poly_commit.random_setup(&mut rand::thread_rng());
        let com_f = poly_commit.commit(&f).unwrap();
        let (com_q, qr, pi_q, fr, pi_f) = gen_proof(w, &f, com_f, degree_w, &poly_commit).unwrap();

        let res = verify(w, com_f, com_q, qr, pi_q, fr, pi_f, degree_w, &poly_commit);
        assert!(res)
    }

    #[test]
    pub fn zero_test_soundness() {
        let w = Fr::from(5u8);
        let degree_w = 3;
        let f = Polynomial::new(vec![
            Fr::from(-375),
            Fr::from(340),
            Fr::from(62),
            Fr::from(-28),
            Fr::from(2),
        ]);
        let mut poly_commit =
            KZG::<Bls12_381>::new(G1Projective::generator(), G2Projective::generator(), 4);
        poly_commit.random_setup(&mut rand::thread_rng());
        let com_f = poly_commit.commit(&f).unwrap();
        let (com_q, qr, pi_q, fr, pi_f) = gen_proof(w, &f, com_f, degree_w, &poly_commit).unwrap();

        let res = verify(w, com_f, com_q, qr, pi_q, fr, pi_f, degree_w, &poly_commit);
        assert!(!res);
    }
}
