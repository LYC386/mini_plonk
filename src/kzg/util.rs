use ark_ff::{Field, PrimeField};

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

#[cfg(test)]

mod test {
    use super::*;
    use ark_bls12_381::Fq;

    #[test]
    fn test_add() {
        //m = 2x^2 + 3x + 4
        let m = vec![Fq::from(4), Fq::from(3), Fq::from(2)];
        // n = x+2
        let n = vec![Fq::from(2), Fq::from(1)];

        let r = add(&m, &n);
        assert_eq!(r, vec![Fq::from(6), Fq::from(4), Fq::from(2)]);
    }

    #[test]
    fn test_mul() {
        //m = 2x^2 + 3x + 4
        let m = vec![Fq::from(4), Fq::from(3), Fq::from(2)];
        // n = x+2
        let n = vec![Fq::from(2), Fq::from(1)];

        let r = mul(&m, &n);
        assert_eq!(r, vec![Fq::from(8), Fq::from(10), Fq::from(7), Fq::from(2)]);
    }

    #[test]
    fn test_eval() {
        //f = 2x^2 + 3x + 4
        let f = vec![Fq::from(4), Fq::from(3), Fq::from(2), Fq::from(0)];

        let r = eval(&f, Fq::from(2));
        assert_eq!(r, Fq::from(18));
    }

    #[test]
    fn test_div() {
        //m = 2x^2 + 3x + 4
        let n = vec![Fq::from(4), Fq::from(3), Fq::from(2), Fq::from(0)];
        // d = x+2
        let d = vec![Fq::from(2), Fq::from(1), Fq::from(0)];

        let (q, r) = div(&n, &d).unwrap();
        assert_eq!(q, vec![Fq::from(-1), Fq::from(2)]);
        assert_eq!(r, vec![Fq::from(6)]);
    }

    #[test]
    fn test_div2() {
        //m = 2x+4
        let n = vec![Fq::from(4), Fq::from(2), Fq::from(0)];
        // d = x+2
        let d = vec![Fq::from(2), Fq::from(1), Fq::from(0)];

        let (q, r) = div(&n, &d).unwrap();
        assert_eq!(q, vec![Fq::from(2)]);
        assert_eq!(r, vec![Fq::from(0)]);
    }

    #[test]
    fn test_div3() {
        // n = x+2
        let n = vec![Fq::from(2), Fq::from(1), Fq::from(0)];
        //d = 2x^2 + 3x + 4
        let d = vec![Fq::from(4), Fq::from(3), Fq::from(2), Fq::from(0)];

        let (q, r) = div(&n, &d).unwrap();
        assert_eq!(q, vec![Fq::from(0)]);
        assert_eq!(r, vec![Fq::from(2), Fq::from(1)]);
    }
}
