use ark_ff::Field;
use num_traits::{One, Zero};
use std::{
    ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub},
    slice::SliceIndex,
    vec,
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

impl<E, Idx> IndexMut<Idx> for Polynomial<E>
where
    Idx: SliceIndex<[E], Output = E>,
{
    fn index_mut(&mut self, index: Idx) -> &mut Self::Output {
        &mut self.data[index]
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

impl<E: Zero + Clone> Add for &Polynomial<E> {
    type Output = Polynomial<E>;

    fn add(self, rhs: &Polynomial<E>) -> Self::Output {
        let mut result = vec![E::zero(); std::cmp::max(self.data.len(), rhs.data.len())];
        for i in 0..result.len() {
            result[i] = self.data.get(i).unwrap_or(&E::zero()).clone()
                + rhs.data.get(i).unwrap_or(&E::zero()).clone()
        }
        Polynomial { data: result }
    }
}

impl<E: Zero + Clone + Sub<Output = E>> Sub for Polynomial<E> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut result = vec![E::zero(); std::cmp::max(self.data.len(), rhs.data.len())];
        for i in 0..result.len() {
            result[i] = self.data.get(i).unwrap_or(&E::zero()).clone()
                - rhs.data.get(i).unwrap_or(&E::zero()).clone()
        }
        Polynomial { data: result }
    }
}

impl<E: Zero + Clone + Sub<Output = E>> Sub for &Polynomial<E> {
    type Output = Polynomial<E>;

    fn sub(self, rhs: &Polynomial<E>) -> Self::Output {
        let mut result = vec![E::zero(); std::cmp::max(self.data.len(), rhs.data.len())];
        for i in 0..result.len() {
            result[i] = self.data.get(i).unwrap_or(&E::zero()).clone()
                - rhs.data.get(i).unwrap_or(&E::zero()).clone()
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

impl<E: Zero + Clone + Mul<Output = E>> Mul for &Polynomial<E> {
    type Output = Polynomial<E>;

    fn mul(self, rhs: &Polynomial<E>) -> Self::Output {
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

impl<E> Div for &Polynomial<E>
where
    E: Zero
        + Clone
        + Div<Output = E>
        + Mul<Output = E>
        + Sub<Output = E>
        + std::cmp::PartialEq
        + std::fmt::Debug,
{
    type Output = Result<(Polynomial<E>, Polynomial<E>), String>;

    fn div(self, rhs: &Polynomial<E>) -> Self::Output {
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
            if *r.last().unwrap() == E::zero() && r.len() != 1 {
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

impl<E> Polynomial<E> {
    pub fn degree(&self) -> usize {
        self.data.len() - 1
    }

    pub fn iter(&self) -> std::slice::Iter<'_, E> {
        self.data.iter()
    }
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
    fn test_sub() {
        //m = 2x^2 + 3x + 4
        let m = Polynomial::new(vec![Fq::from(4), Fq::from(3), Fq::from(2)]);
        // n = x+2
        let n = Polynomial::new(vec![Fq::from(2), Fq::from(1)]);

        let r = n - m;
        assert_eq!(
            r,
            Polynomial::new(vec![Fq::from(-2), Fq::from(-2), Fq::from(-2)])
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

        let (q, r) = (&n / &d).unwrap();
        assert_eq!(q, Polynomial::new(vec![Fq::from(0)]));
        assert_eq!(r, Polynomial::new(vec![Fq::from(2), Fq::from(1)]));
    }

    #[test]
    fn test_div4() {
        // n = 2x+3
        let n = Polynomial::new(vec![7, 0, 0, 0, 0, 0, 0, 0, 1]);
        //d = 2x^4 + 3x^3 + 0x^2 + x + 4
        let d = Polynomial::new(vec![7, 8, 0, 7, 0, 0, 0, 0, 4, 3, 0, 2]);

        let (q, r) = (&d / &n).unwrap();
        // assert_eq!(q, Polynomial::new(vec![0]));
        // assert_eq!(r, Polynomial::new(vec![Fq::from(2), Fq::from(1)]));
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
