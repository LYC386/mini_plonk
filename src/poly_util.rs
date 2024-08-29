use ark_ff::Field;
use num_traits::{One, Pow, Zero};
use std::{
    ops::{Add, Div, Index, Mul, Neg, Sub},
    slice::SliceIndex,
};

#[derive(Debug, Clone)]
pub struct Polynomial<E> {
    data: Vec<E>,
}

impl<E> Polynomial<E> {
    pub fn new(coeff: Vec<E>) -> Self {
        Polynomial { data: coeff }
    }
}

impl<E: Field> Polynomial<E> {
    pub fn eval(&self, v: E) -> E {
        let mut res = E::zero();
        for (f_i, i) in self.data.iter().zip(0u64..) {
            res = res + f_i.clone() * v.pow(&[i]);
        }
        res
    }
}

impl<E: PartialEq> PartialEq for Polynomial<E> {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

impl<E> Polynomial<E>
where
    E: Clone + Zero + One + Neg<Output = E> + Sub<Output = E> + Div<Output = E>,
{
    pub fn interpolate(points: &[(E, E)]) -> Self {
        let mut result = Polynomial::new(vec![E::zero()]);
        for i in 0..points.len() {
            let mut l_i = Polynomial::new(vec![E::one()]);
            for j in 0..points.len() {
                if i == j {
                    continue;
                }
                let d = points[i].0.clone() - points[j].0.clone();
                let l_i_j = Polynomial::new(vec![-points[j].0.clone() / d.clone(), E::one() / d]);
                l_i = l_i * l_i_j;
            }
            // result = add(&result, &mul(&[points[i].1], &l_i));
            result = result + (Polynomial::new(vec![points[i].1.clone()]) * l_i)
        }
        return result;
    }
}

impl<E, Idx> Index<Idx> for Polynomial<E>
where
    Idx: SliceIndex<[E], Output = E>,
{
    type Output = E;

    fn index(&self, index: Idx) -> &Self::Output {
        &self.data[index]
    }
}

impl<E: Zero + Clone> Add for Polynomial<E> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut result = vec![E::zero(); std::cmp::max(self.data.len(), rhs.data.len())];
        for i in 0..result.len() {
            result[i] = self.data.get(i).unwrap_or(&E::zero()).clone()
                + rhs.data.get(i).unwrap_or(&E::zero()).clone()
        }
        Polynomial { data: result }
    }
}

impl<E: Zero + Clone + Mul<Output = E>> Mul for Polynomial<E> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = vec![E::zero(); self.data.len() + rhs.data.len() - 1];
        for (pos_i, i) in self.data.iter().enumerate() {
            for (pos_j, j) in rhs.data.iter().enumerate() {
                result[pos_i + pos_j] = result[pos_i + pos_j].clone() + i.clone() * j.clone();
            }
        }
        Polynomial { data: result }
    }
}

//polynomial long division, returns (quotient, remainder)
impl<E> Div for Polynomial<E>
where
    E: Zero + Clone + Div<Output = E> + Mul<Output = E> + Sub<Output = E> + std::cmp::PartialEq,
{
    type Output = Result<(Self, Self), String>;

    fn div(self, rhs: Self) -> Self::Output {
        // check empty list
        if self.data.len() == 0 || rhs.data.len() == 0 {
            return Err("empty list".into());
        }

        // check divided by zero
        for i in &rhs.data {
            if *i != E::zero() {
                break;
            }
            return Err("divide by zero".into());
        }

        let mut r: Vec<E> = self.data.clone();
        let mut d: Vec<E> = rhs.data.clone();

        // clean leading zeros
        while *r.last().unwrap() == E::zero() {
            if r.len() != 1 {
                r.pop();
            }
        }
        while *d.last().unwrap() == E::zero() {
            // already checked non-zero
            d.pop();
        }

        if r.len() < d.len() {
            return Ok((
                Polynomial {
                    data: vec![E::zero()],
                },
                Polynomial { data: r },
            ));
        }

        // calculate quotient degree
        let mut q = vec![E::zero(); r.len() - d.len() + 1];
        let mut q_pos = q.len();

        while r.len() >= d.len() {
            q_pos -= 1;
            // t = lead(r)/lead(d)
            let t = r.last().unwrap().clone() / d.last().unwrap().clone();
            // q = q + t
            q[q_pos] = t.clone();

            //r = r - t * d
            for (i, di) in d.iter().enumerate() {
                r[i + q_pos] = r[i + q_pos].clone() - t.clone() * di.clone();
            }

            // clean leading zeros of r
            while *r.last().unwrap() == E::zero() && r.len() != 1 {
                r.pop();
            }
            // if r=0, return
            if r.len() == 1 && r[0] == E::zero() {
                break;
            }
        }
        Ok((Polynomial { data: q }, Polynomial { data: r }))
    }
}

impl<E: Zero + std::cmp::PartialEq> Polynomial<E> {
    pub fn is_zero(&self) -> bool {
        for i in &self.data {
            if *i != E::zero() {
                return false;
            }
        }
        true
    }
}
//polynomial addition
pub fn add<E: Field>(n: &[E], m: &[E]) -> Vec<E> {
    let mut result = vec![E::ZERO; std::cmp::max(n.len(), m.len())];
    for i in 0..result.len() {
        result[i] = *n.get(i).unwrap_or(&E::ZERO) + m.get(i).unwrap_or(&E::ZERO)
    }
    result
}

//polynomial multiplication
pub fn mul<E: Field>(n: &[E], m: &[E]) -> Vec<E> {
    let mut result = vec![E::ZERO; n.len() + m.len() - 1];
    for (pos_i, i) in n.iter().enumerate() {
        for (pos_j, j) in m.iter().enumerate() {
            result[pos_i + pos_j] += *i * j;
        }
    }
    result
}

// evaluate a polynomial f on v
pub fn eval<E: Field>(f: &[E], v: E) -> E {
    let mut res = E::ZERO;
    for (f_i, i) in f.iter().zip(0u64..) {
        res += *f_i * v.pow(&[i]);
    }
    res
}

//polynomial long division
// n/d returns (quotient, remainder)
pub fn div<E: Field>(n: &[E], d: &[E]) -> Result<(Vec<E>, Vec<E>), String> {
    // check empty list
    if d.len() == 0 || n.len() == 0 {
        return Err("empty list".into());
    }

    // check divided by zero
    for i in d {
        if *i != E::ZERO {
            break;
        }
        return Err("divide by zero".into());
    }

    let mut r: Vec<E> = n.to_vec();
    let mut d: Vec<E> = d.to_vec();

    // clean leading zeros
    while *r.last().unwrap() == E::ZERO {
        if r.len() != 1 {
            r.pop();
        }
    }
    while *d.last().unwrap() == E::ZERO {
        // already checked non-zero
        d.pop();
    }

    if r.len() < d.len() {
        return Ok((vec![E::ZERO], r));
    }

    // calculate quotient degree
    let mut q = vec![E::ZERO; r.len() - d.len() + 1];
    let mut q_pos = q.len();

    while r.len() >= d.len() {
        q_pos -= 1;
        // t = lead(r)/lead(d)
        let t = *r.last().unwrap() / d.last().unwrap();
        // q = q + t
        q[q_pos] = t;

        //r = r - t * d
        for (i, di) in d.iter().enumerate() {
            r[i + q_pos] -= t * di;
        }

        // clean leading zeros of r
        while *r.last().unwrap() == E::ZERO && r.len() != 1 {
            r.pop();
        }
        // if r=0, return
        if r.len() == 1 && r[0] == E::ZERO {
            break;
        }
    }
    Ok((q, r))
}

pub fn is_zero<E: Field>(f: &[E]) -> bool {
    for i in f {
        if *i != E::ZERO {
            return false;
        }
    }
    true
}

// interpolate a function with input points (x,y) with Lagrange interpolation
pub fn interpolate<E: Field>(points: &[(E, E)]) -> Vec<E> {
    let mut result = vec![E::ZERO];
    for i in 0..points.len() {
        let mut l_i = vec![E::ONE];
        for j in 0..points.len() {
            if i == j {
                continue;
            }
            let d = points[i].0 - points[j].0;
            let l_i_j = &[-points[j].0 / d, E::ONE / d];
            l_i = mul(&l_i, l_i_j);
        }
        result = add(&result, &mul(&[points[i].1], &l_i));
    }
    return result;
}

#[cfg(test)]

mod test {
    use super::*;
    use ark_bls12_381::Fq;

    #[test]
    fn test_add() {
        //m = 2x^2 + 3x + 4
        let m = Polynomial::new(vec![Fq::from(4), Fq::from(3), Fq::from(2)]);
        // n = x+2
        let n = Polynomial::new(vec![Fq::from(2), Fq::from(1)]);

        let r = m + n;
        assert_eq!(
            r,
            Polynomial::new(vec![Fq::from(6), Fq::from(4), Fq::from(2)])
        );
    }

    #[test]
    fn test_mul() {
        //m = 2x^2 + 3x + 4
        let m = Polynomial::new(vec![Fq::from(4), Fq::from(3), Fq::from(2)]);
        // n = x+2
        let n = Polynomial::new(vec![Fq::from(2), Fq::from(1)]);

        let r = m * n;
        assert_eq!(
            r,
            Polynomial::new(vec![Fq::from(8), Fq::from(10), Fq::from(7), Fq::from(2)])
        );

        let k = Polynomial::new(vec![Fq::from(2)]);
        let r2 = r * k;
        assert_eq!(
            r2,
            Polynomial::new(vec![Fq::from(16), Fq::from(20), Fq::from(14), Fq::from(4)])
        );
    }

    #[test]
    fn test_eval() {
        //f = 2x^2 + 3x + 4
        let f = Polynomial::new(vec![Fq::from(4), Fq::from(3), Fq::from(2), Fq::from(0)]);

        let r = f.eval(Fq::from(2));
        assert_eq!(r, Fq::from(18));
    }

    #[test]
    fn test_div() {
        //m = 2x^2 + 3x + 4
        let n = Polynomial::new(vec![Fq::from(4), Fq::from(3), Fq::from(2), Fq::from(0)]);
        // d = x+2
        let d = Polynomial::new(vec![Fq::from(2), Fq::from(1), Fq::from(0)]);

        let (q, r) = (n / d).unwrap();
        assert_eq!(q, Polynomial::new(vec![Fq::from(-1), Fq::from(2)]));
        assert_eq!(r, Polynomial::new(vec![Fq::from(6)]));
    }

    #[test]
    fn test_div2() {
        //m = 2x+4
        let n = Polynomial::new(vec![Fq::from(4), Fq::from(2), Fq::from(0)]);
        // d = x+2
        let d = Polynomial::new(vec![Fq::from(2), Fq::from(1), Fq::from(0)]);

        let (q, r) = (n / d).unwrap();
        assert_eq!(q, Polynomial::new(vec![Fq::from(2)]));
        assert_eq!(r, Polynomial::new(vec![Fq::from(0)]));
    }

    #[test]
    fn test_div3() {
        // n = x+2
        let n = Polynomial::new(vec![Fq::from(2), Fq::from(1), Fq::from(0)]);
        //d = 2x^2 + 3x + 4
        let d = Polynomial::new(vec![Fq::from(4), Fq::from(3), Fq::from(2), Fq::from(0)]);

        let (q, r) = (n / d).unwrap();
        assert_eq!(q, Polynomial::new(vec![Fq::from(0)]));
        assert_eq!(r, Polynomial::new(vec![Fq::from(2), Fq::from(1)]));
    }

    #[test]

    fn test_interpolate() {
        //f = 2x^2 + 3x + 4
        let points = &[
            (Fq::from(0), Fq::from(4)),
            (Fq::from(1), Fq::from(9)),
            (Fq::from(3), Fq::from(31)),
        ];
        let res = Polynomial::interpolate(points);
        assert_eq!(
            res,
            Polynomial::new(vec![Fq::from(4), Fq::from(3), Fq::from(2)])
        );
    }
}
