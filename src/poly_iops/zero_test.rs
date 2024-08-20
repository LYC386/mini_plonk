use crate::{kzg::KZG, poly_util};
use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use crypto_hash::{Algorithm::SHA256, Hasher};
use std::{io::Write, ops::Neg};

// generate a zero proof, returns com_q, q(r), pi_q, f(r), pi_f
pub fn gen_proof<E: Pairing>(
    w: E::ScalarField,
    f: &[E::ScalarField],
    com_f: E::G1,
    degree_w: usize,
    poly_commit: KZG<E>,
) -> Result<(E::G1, E::ScalarField, E::G1, E::ScalarField, E::G1), String> {
    // vanishing polynimial of W = P{1, w, w^2, ..., w^k-1}
    let mut z_w = vec![E::ScalarField::ZERO; degree_w];
    z_w[0] = E::ScalarField::from(1u32).neg();
    z_w[degree_w - 1] = E::ScalarField::from(1u32);

    // calaulate q(x)
    let (q, r) = poly_util::div(f, &z_w)?;
    if !poly_util::is_zero(&r) {
        return Err("remainder is not zero".into());
    }
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
    poly_commit: KZG<E>,
) -> bool {
    // vanishing polynimial of W = P{1, w, w^2, ..., w^k-1}
    let mut z_w = vec![E::ScalarField::ZERO; degree_w];
    z_w[0] = E::ScalarField::from(1u32).neg();
    z_w[degree_w - 1] = E::ScalarField::from(1u32);

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
    let z_wr = poly_util::eval(&z_w, r);
    if fr != qr * z_wr {
        return false;
    }
    true
}
