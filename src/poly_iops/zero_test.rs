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
    // vanishing polynimial of W = {1, w, w^2, ..., w^k-1} where w is kth root of unity
    let mut z_w = vec![E::ScalarField::ZERO; degree_w + 1];
    z_w[0] = -E::ScalarField::from(1u32);
    z_w[degree_w] = E::ScalarField::from(1u32);
    let z_w = Polynomial::new(z_w);

    // calaulate q(x)
    let (q, _) = (f / &z_w)?;
    let com_q = poly_commit.commit(&q)?;

    // generate r
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    hasher.write_all(b"ZERO_TEST").unwrap();
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
    // vanishing polynimial of W = {1, w, w^2, ..., w^k-1} where w is kth root of unity
    let mut z_w = vec![E::ScalarField::ZERO; degree_w + 1];
    z_w[0] = -E::ScalarField::from(1u32);
    z_w[degree_w] = E::ScalarField::from(1u32);
    let z_w = Polynomial::new(z_w);

    // calculate coin r
    let mut hasher = Hasher::new(SHA256);
    w.serialize_uncompressed(&mut hasher).unwrap();
    com_f.serialize_uncompressed(&mut hasher).unwrap();
    hasher.write_all(&degree_w.to_be_bytes()).unwrap();
    hasher.write_all(b"ZERO_TEST").unwrap();
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
    use ark_ff::PrimeField;
    use ark_std::One;
    use num_bigint::BigUint;

    #[test]
    pub fn zero_test_happy_path() {
        // create w as an 8th root of unity
        let modulus: BigUint = Fr::MODULUS.into();
        let one: num_bigint::BigUint = Fr::one().into();
        let ord = &modulus - one;
        let fr_gen = Fr::from(7);
        let subgroup_ord = Fr::from(8);
        let w = fr_gen.pow((ord / BigUint::from(subgroup_ord)).to_u64_digits());
        let degree_w = 8;

        // f = (x^8-1)(2x^3+3x+4)
        let f = Polynomial::new(vec![
            Fr::from(-1),
            Fr::ZERO,
            Fr::ZERO,
            Fr::ZERO,
            Fr::ZERO,
            Fr::ZERO,
            Fr::ZERO,
            Fr::ZERO,
            Fr::from(1),
        ]);
        let f = f * Polynomial::new(vec![Fr::from(4), Fr::from(3), Fr::from(0), Fr::from(2)]);
        let mut poly_commit = KZG::<Bls12_381>::new(
            G1Projective::generator(),
            G2Projective::generator(),
            f.degree(),
        );
        poly_commit.random_setup(&mut rand::thread_rng());
        let com_f = poly_commit.commit(&f).unwrap();
        let (com_q, qr, pi_q, fr, pi_f) = gen_proof(w, &f, com_f, degree_w, &poly_commit).unwrap();

        let res = verify(w, com_f, com_q, qr, pi_q, fr, pi_f, degree_w, &poly_commit);
        assert!(res)
    }

    #[test]
    pub fn zero_test_soundness() {
        // create w as an 8th root of unity
        let modulus: BigUint = Fr::MODULUS.into();
        let one: num_bigint::BigUint = Fr::one().into();
        let ord = &modulus - one;
        let fr_gen = Fr::from(7);
        let subgroup_ord = Fr::from(8);
        let w = fr_gen.pow((ord / BigUint::from(subgroup_ord)).to_u64_digits());
        let degree_w = 8;

        // f = (x^8-1)(2x^3+3x+4)
        let f = Polynomial::new(vec![
            Fr::from(-1),
            Fr::ZERO,
            Fr::ONE,
            Fr::ZERO,
            Fr::ZERO,
            Fr::ZERO,
            Fr::ZERO,
            Fr::ZERO,
            Fr::from(1),
        ]);
        let f = f * Polynomial::new(vec![Fr::from(4), Fr::from(3), Fr::from(0), Fr::from(2)]);
        let mut poly_commit = KZG::<Bls12_381>::new(
            G1Projective::generator(),
            G2Projective::generator(),
            f.degree(),
        );
        poly_commit.random_setup(&mut rand::thread_rng());
        let com_f = poly_commit.commit(&f).unwrap();
        let (com_q, qr, pi_q, fr, pi_f) = gen_proof(w, &f, com_f, degree_w, &poly_commit).unwrap();

        let res = verify(w, com_f, com_q, qr, pi_q, fr, pi_f, degree_w, &poly_commit);
        assert!(!res)
    }
}
