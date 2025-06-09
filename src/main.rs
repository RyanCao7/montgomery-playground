//! Ryan playground implementation of Montgomery arithmetic.
//!
//! Implicitly, we are always using R = 2^{16} as an example here.

use num_primes::Generator;
use rand::RngCore;

use crate::utils::{MOD_U32_BIT_MASK, compute_inverse, div_r_u16, mod_r_u16, wide_mul_u16};
pub mod utils;

/// Precomputed set of constants for each modulus which are used within
/// Montgomery multiplication + conversion.
///
/// Note: `n_prime` = -(n^{-1}) \mod R, such that N N' = -1 \mod R
struct MontgomeryParametersU16 {
    modulus: u16,
    r_squared_mod_n: u16,
    n_prime: u32,
}

impl MontgomeryParametersU16 {
    pub fn new(modulus: u16) -> Self {
        let r = (1 as u32) << 16;
        let r_mod_n = r % (modulus as u32);
        let r_squared_mod_n = ((r_mod_n * r_mod_n) % (modulus as u32))
            .try_into()
            .expect("This should work");

        Self {
            modulus,
            r_squared_mod_n,
            n_prime: r - compute_inverse(modulus as u32, r),
        }
    }
}

/// Brings a value 0 \leq a < N into its "Montgomery form", aR \mod N.
/// We assume here that R = 2^{32} and that N < R.
///
/// The trick is to simply perform one wide mul by R^2 \mod N and then perform
/// a single reduction step, since now we are in the residue class aR^2 \mod N.
fn convert_to_montgomery_form_u16(a: u16, montgomery_parameters: &MontgomeryParametersU16) -> u16 {
    let a_times_r_squared = wide_mul_u16(a, montgomery_parameters.r_squared_mod_n);

    // This should get us to just aR \mod N
    montgomery_reduction_u16(a_times_r_squared, montgomery_parameters)
}

/// Computes TR^{-1} \mod N. Note that R = 2^{32} in this case.
///
/// ## Formula:
/// (T + ((T \mod R N') \mod R) N) / R
///
/// ## Assumptions
/// * 0 \leq T < RN
fn montgomery_reduction_u16(t_mont: u32, montgomery_parameters: &MontgomeryParametersU16) -> u16 {
    debug_assert!(t_mont < ((1 as u32) << 16) * (montgomery_parameters.modulus as u32));
    let t_mod_r = mod_r_u16(t_mont);
    // Note: \alpha = ((T \mod R N') \mod R).
    let alpha = mod_r_u16((t_mod_r as u32) * montgomery_parameters.n_prime);
    // We compute (T + \alpha N) / R
    let numerator = t_mont + wide_mul_u16(alpha, montgomery_parameters.modulus);
    // The lower 16 bits should all be zero, i.e. divisible by 2^{16}
    debug_assert!(numerator & MOD_U32_BIT_MASK == 0);
    div_r_u16(numerator) // Finally, right-shift by R
}

fn main() {
    let modulus_digits = Generator::safe_prime(16).to_u32_digits();
    debug_assert_eq!(modulus_digits.len(), 1);
    let modulus = modulus_digits[0] as u16;
    dumb_prime_checker(modulus);

    let mut rng = rand::rng();
    let a = (rng.next_u32() >> 16) as u16;
    let b = (rng.next_u32() >> 16) as u16;
    let raw_ab_product = wide_mul_u16(a, b); // Full width multiplication
    let ab_mod_modulus: u16 = (raw_ab_product % (modulus as u32))
        .try_into()
        .expect("The residue class should fit in a u16");

    let params = MontgomeryParametersU16::new(modulus);
    let a_mont_form = convert_to_montgomery_form_u16(a, &params);
    let b_mont_form = convert_to_montgomery_form_u16(b, &params);

    let a_times_b_mont_form =
        montgomery_reduction_u16(wide_mul_u16(a_mont_form, b_mont_form), &params);
    let ab_mod_modulus_via_montgomery =
        montgomery_reduction_u16(a_times_b_mont_form as u32, &params);

    assert_eq!(ab_mod_modulus, ab_mod_modulus_via_montgomery);
}

/// Simple checker that ensures no divisibility up to sqrt(n)
fn dumb_prime_checker(n: u16) {
    let upper_bound = ((n as f64).sqrt()).ceil() as u16;
    (2..(upper_bound + 1)).for_each(|val| {
        if n % val == 0 {
            dbg!(n);
            dbg!(val);
        }
        assert_ne!(n % val, 0);
    });
}

#[cfg(test)]
pub mod tests {
    use std::time::Instant;

    use num_primes::Generator;
    use rand::RngCore;

    use crate::{
        MontgomeryParametersU16, convert_to_montgomery_form_u16, dumb_prime_checker,
        montgomery_reduction_u16,
        utils::{wide_mul_u16, wide_mul_u32},
    };

    /// Okay well this is a complete failure but oh well LOL
    #[test]
    pub fn basic_mult_vs_mod_timing_test() {
        let mut test_rng = rand::rng();
        let rand_a: Vec<u32> = (0..1000000).map(|_| test_rng.next_u32()).collect();
        let rand_b: Vec<u32> = (0..1000000).map(|_| test_rng.next_u32()).collect();
        // let shift_b: Vec<u8> = (0..1000000).map(|_| test_rng.random_range(0..30)).collect();

        let mult_start = Instant::now();
        rand_a
            .into_iter()
            .zip(rand_b.into_iter())
            .for_each(|(a, b)| {
                let _ = wide_mul_u32(a, b);
            });
        let mult_timer = mult_start.elapsed();

        let other_rand_a: Vec<u32> = (0..1000000).map(|_| test_rng.next_u32()).collect();
        let other_rand_b: Vec<u32> = (0..1000000).map(|_| test_rng.next_u32()).collect();
        let mod_start = Instant::now();
        other_rand_a
            .into_iter()
            .zip(other_rand_b.into_iter())
            .for_each(|(a, b)| {
                let _ = a % b;
            });
        let mod_timer = mod_start.elapsed();

        dbg!(mult_timer);
        dbg!(mod_timer);
    }

    /// The answer is probably "no" lol
    #[test]
    fn does_montgomery_mul_actually_speed_things_up() {
        let modulus_digits = Generator::safe_prime(16).to_u32_digits();
        assert_eq!(modulus_digits.len(), 1);
        let modulus = modulus_digits[0] as u16;
        dumb_prime_checker(modulus);

        let mut test_rng = rand::rng();
        let rand_a: Vec<u16> = (0..10000000)
            .map(|_| (test_rng.next_u32() >> 16) as u16)
            .collect();
        let rand_a_copy = rand_a.clone();

        // Let's do the naive version first
        let naive_start = Instant::now();
        let naive_result = rand_a
            .into_iter()
            .reduce(|acc, a| {
                let raw_wide_product = wide_mul_u16(acc, a);
                (raw_wide_product % (modulus as u32))
                    .try_into()
                    .expect("The residue class should fit in a u16")
            })
            .expect("This should not fail. Lol.");
        let naive_time = naive_start.elapsed();

        let params = MontgomeryParametersU16::new(modulus);

        let montgomery_start = Instant::now();
        let first_val = rand_a_copy[0];
        let first_val_montgomery_form = convert_to_montgomery_form_u16(first_val, &params);
        let montgomery_result =
            rand_a_copy
                .into_iter()
                .fold(first_val_montgomery_form, |acc, a| {
                    // First, we convert `a` into montgomery form
                    let a_mont_form = convert_to_montgomery_form_u16(a, &params);
                    // Then, we perform a wide multiply
                    let wide_product = wide_mul_u16(acc, a_mont_form);
                    // Finally, a REDC operation
                    montgomery_reduction_u16(wide_product, &params)
                });
        // One final reduction to take things fully out of Montgomery form
        let final_montgomery_result = montgomery_reduction_u16(montgomery_result as u32, &params);
        let montgomery_time = montgomery_start.elapsed();

        assert_eq!(naive_result, final_montgomery_result);

        dbg!(naive_time);
        dbg!(montgomery_time);
    }
}
