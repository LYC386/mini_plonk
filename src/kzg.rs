mod util;

use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_std::UniformRand;
use ark_std::Zero;
use rand::Rng;
use std::ops::Mul;
use std::ops::Neg;
pub struct KZG<E: Pairing> {
    pub g1: E::G1,
    pub g2: E::G2,
    pub gp1: Vec<E::G1>,
    pub gp2: Vec<E::G2>,
    pub degree: usize,
}

impl<E: Pairing> KZG<E> {
    pub fn new(g1: E::G1, g2: E::G2, degree: usize) -> Self {
        KZG {
            g1,
            g2,
            gp1: vec![],
            gp2: vec![],
            degree,
        }
    }

    // generate a random t and
    // create gp1 & gp2 as [g1, t*g1, t^2*g1, ..., t^d*g1]
    // and [g2, t*g2, t^2*g2, ..., t^d*g2]
    pub fn random_setup<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        let tau = E::ScalarField::rand(rng);
        for i in 0..self.degree + 1 {
            self.gp1.push(self.g1.mul(tau.pow(&[i as u64])));
            self.gp2.push(self.g2.mul(tau.pow(&[i as u64])));
        }
    }

    // f_0 * gp1_0 +...+f_d * gp1_d
    // calculate commitment of f (in coefficient repre)
    pub fn commit(&self, f: &[E::ScalarField]) -> Result<E::G1, String> {
        let mut com_f = E::G1::zero();
        if f.len() > self.degree {
            return Err("degree of f exceeds maximum degree".into());
        }
        for (pi, fi) in self.gp1.iter().zip(f.iter()) {
            com_f = com_f + pi.mul(fi);
        }
        Ok(com_f)
    }

    // evaluate f(u) and return (f(u), pi)
    pub fn eval(
        &self,
        f: &[E::ScalarField],
        u: E::ScalarField,
    ) -> Result<(E::ScalarField, E::G1), String> {
        let mut f = f.to_vec();
        let f_u = util::eval(&f, u);
        // f(x) = f(x) - f(u)
        f[0] = f[0] - f_u;

        // d(x) = (x-u)
        let dx = vec![u.neg(), E::ScalarField::from(1u64)];

        // q(x) = f(x)-f(u) / d(x)
        let (qx, r) = util::div(&f, &dx)?;
        if !util::is_zero(&r) {
            return Err("r(x) is not zero".into());
        }

        // pi = q(t) * g
        let pi = self.commit(&qx)?;
        return Ok((f_u, pi));
    }

    pub fn verify(&self, com_f: E::G1, u: E::ScalarField, v: E::ScalarField, pi: E::G1) -> bool {
        let tg2_m_ug2: E::G2 = self.gp2[1] - self.g2 * u;
        let comf_m_vg1: E::G1 = com_f - self.g1 * v;

        E::pairing(pi, tg2_m_ug2) == E::pairing(comf_m_vg1, self.g2)
    }
}

#[cfg(test)]

mod test {
    use super::*;
    use ark_bls12_381::{Bls12_381, Fr, G1Projective, G2Projective};
    use ark_ec::Group;

    pub fn test_setup() -> KZG<Bls12_381> {
        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();
        let mut kzg = KZG::<Bls12_381>::new(g1, g2, 11);
        let mut rng = rand::thread_rng();
        kzg.random_setup(&mut rng);

        kzg
    }

    #[test]
    pub fn test_e2e_happy_path() {
        let kzg = test_setup();
        //f = 2x^2 + 3x + 4
        let f = vec![Fr::from(4), Fr::from(3), Fr::from(2), Fr::from(0)];
        let com_f = kzg.commit(&f).unwrap();
        let u = Fr::from(2);
        let (v, pi) = kzg.eval(&f, u).unwrap();
        assert!(v == Fr::from(18));
        let r = kzg.verify(com_f, u, v, pi);
        assert!(r);
    }

    #[test]
    pub fn test_wrong_pi() {
        let kzg = test_setup();
        //f = 2x^2 + 3x + 4
        let f = vec![Fr::from(4), Fr::from(3), Fr::from(2), Fr::from(0)];
        let f_prime = vec![Fr::from(4), Fr::from(3), Fr::from(1)];
        let com_f = kzg.commit(&f_prime).unwrap();
        let u = Fr::from(2);
        let (v, pi) = kzg.eval(&f, u).unwrap();
        assert!(v == Fr::from(18));
        let r = kzg.verify(com_f, u, v, pi);
        assert!(!r);
    }
}
