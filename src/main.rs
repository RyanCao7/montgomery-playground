//! Ryan playground implementation of Montgomery arithmetic.
//!
//! Implicitly, we are always using R = 2^{16} as an example here.

use num_primes::Generator;
use rand::RngCore;

use crate::{
    montgomery::{
        MontgomeryParametersU16, convert_to_montgomery_form_u16, montgomery_reduction_u16,
    },
    utils::{dumb_prime_checker, wide_mul_u16},
};
pub mod barrett;
pub mod montgomery;
pub mod utils;

fn main() {
    // --- Montgomery multiplication example (single mul) ---
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

#[cfg(test)]
pub mod tests {
    use std::time::Instant;

    use num_primes::Generator;
    use rand::RngCore;

    use crate::{
        MontgomeryParametersU16,
        barrett::{BarrettParamsU16, barrett_reduction_u16},
        convert_to_montgomery_form_u16, montgomery_reduction_u16,
        utils::{dumb_prime_checker, wide_mul_u16, wide_mul_u32},
    };

    /// running 1 test
    /// [src/main.rs:147:9] mult_timer = 122.125µs
    /// [src/main.rs:148:9] mod_timer = 491.625µs
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

    /// running 1 test
    /// [src/main.rs:196:9] naive_time = 384.274709ms
    /// [src/main.rs:197:9] montgomery_time = 339.6625ms
    /// test tests::does_montgomery_mul_actually_speed_things_up ... ok
    ///
    /// TODO(ryancao): Do this but in Criterion
    #[test]
    fn does_montgomery_mul_actually_speed_things_up() {
        let modulus_digits = Generator::safe_prime(16).to_u32_digits();
        assert_eq!(modulus_digits.len(), 1);
        let modulus = modulus_digits[0] as u16;
        dumb_prime_checker(modulus);

        let mut test_rng = rand::rng();
        let rand_a: Vec<u16> = (0..100000000)
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

    /// running 1 test
    /// [src/main.rs:242:9] naive_time = 415.469833ms
    /// [src/main.rs:243:9] barrett_time = 337.670083ms
    /// test tests::test_barrett ... ok
    ///
    /// TODO(ryancao): Do this but in Criterion
    #[test]
    fn test_barrett() {
        let modulus_digits = Generator::safe_prime(16).to_u32_digits();
        assert_eq!(modulus_digits.len(), 1);
        let modulus = modulus_digits[0] as u16;
        dumb_prime_checker(modulus);

        let mut test_rng = rand::rng();
        let rand_a: Vec<u16> = (0..100000000)
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

        let barrett_params = BarrettParamsU16::new(modulus);

        let barrett_start = Instant::now();
        let barrett_result = rand_a_copy
            .into_iter()
            .reduce(|acc, a| {
                // First we perform a wide multiply
                let naive_mul_no_reduction = wide_mul_u16(acc, a);
                // Then a Barrett reduction
                barrett_reduction_u16(naive_mul_no_reduction, &barrett_params)
            })
            .expect("This should work!");
        let barrett_time = barrett_start.elapsed();

        assert_eq!(naive_result, barrett_result);

        dbg!(naive_time);
        dbg!(barrett_time);
    }
}
