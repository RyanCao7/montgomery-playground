/// Returns gcd(a, b) using Euclid's algorithm.
pub fn compute_euclid_gcd(a: u32, b: u32) -> u32 {
    if b == 0 {
        return a;
    }
    compute_euclid_gcd(b, a % b)
}

/// ax + by = gcd(a, b)
#[derive(Clone, Debug, Copy)]
pub struct ExtendedEuclidResult {
    a: u32,
    b: u32,
    gcd: u32,
    x: i128,
    y: i128,
}

pub fn extended_euclid_gcd_bezout(a: u32, b: u32) -> ExtendedEuclidResult {
    // Make sure to remove this for release mode
    let euclid_gcd = compute_euclid_gcd(a, b);
    dbg!(euclid_gcd);

    let (a_internal, b_internal) = if a > b {
        (a as i128, b as i128)
    } else {
        (b as i128, a as i128)
    };

    let mut r_i_minus_two: i128 = a_internal;
    let mut r_i_minus_one: i128 = b_internal;
    let mut s_i_minus_two: i128 = 1;
    let mut s_i_minus_one: i128 = 0;
    let mut t_i_minus_two: i128 = 0;
    let mut t_i_minus_one: i128 = 1;

    while r_i_minus_one != 0 {
        let q_i_minus_two = r_i_minus_two / r_i_minus_one;
        // r_i = r_{i - 2} - r_{i - 1} * q_{i - 2}
        let r_i = r_i_minus_two - r_i_minus_one * q_i_minus_two;
        // s_i = (s_{i - 2} - s_{i - 1} * q_{i - 2})
        let s_i = s_i_minus_two - s_i_minus_one * q_i_minus_two;
        // t_i = (t_{i - 2} - t_{i - 1} * q_{i - 2})
        let t_i = t_i_minus_two - t_i_minus_one * q_i_minus_two;

        (r_i_minus_one, r_i_minus_two) = (r_i, r_i_minus_one);
        (s_i_minus_one, s_i_minus_two) = (s_i, s_i_minus_one);
        (t_i_minus_one, t_i_minus_two) = (t_i, t_i_minus_one);
    }

    // Final sanitychecks to ensure we actually found the correct GCD
    // and (x, y) values for Bezout's lemma
    let extended_euclid_gcd = r_i_minus_two as u32;
    let (x, y) = if a > b {
        (s_i_minus_two, t_i_minus_two)
    } else {
        (t_i_minus_two, s_i_minus_two)
    };
    debug_assert_eq!(
        extended_euclid_gcd as i128,
        (a as i128) * x + (b as i128) * y
    );
    debug_assert_eq!(extended_euclid_gcd, euclid_gcd);

    ExtendedEuclidResult {
        a: a as u32,
        b: b as u32,
        gcd: extended_euclid_gcd,
        x,
        y,
    }
}

pub fn wide_mul_u32(a: u32, b: u32) -> u64 {
    (a as u64) * (b as u64)
}

pub fn wide_mul_u16(a: u16, b: u16) -> u32 {
    (a as u32) * (b as u32)
}

/// Computes x^{-1} in the field modulo `n`.
/// Note: Assumes that `x` and `n` are relatively prime!
pub fn compute_inverse(x: u32, n: u32) -> u32 {
    let euclid_result = extended_euclid_gcd_bezout(x, n);
    dbg!(euclid_result);
    dbg!((euclid_result.a as i128) * euclid_result.x + (euclid_result.b as i128) * euclid_result.y);
    let inv = (euclid_result.x % (n as i128)) as u32;
    // let inv = if euclid_result.x < 0 {
    //     ((-euclid_result.x) % (n as i128)) as u32
    // } else {
    //     (euclid_result.x % (n as i128)) as u32
    // };

    debug_assert_eq!(wide_mul_u32(inv, x) % (n as u64), 1);

    inv
}

#[cfg(test)]
pub mod tests {
    use crate::utils::extended_euclid_gcd_bezout;
    use rand::{TryRngCore, rngs::OsRng};

    #[test]
    fn test_extended_euclid() {
        let mut rng = OsRng::default();
        (0..100).for_each(|idx| {
            let rand_a = rng.try_next_u32().expect("Should've been able to generate");
            let rand_b = rng.try_next_u32().expect("Should've been able to generate");
            let _ = extended_euclid_gcd_bezout(rand_a, rand_b);
            dbg!(idx);
        });
    }
}
