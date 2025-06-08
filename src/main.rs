//! Ryan playground implementation of Montgomery arithmetic.
//!
//! Implicitly, we are always using R = 2^{32} as an example here.

use std::convert;

use num_primes::Generator;
use rand::RngCore;

use crate::utils::{compute_inverse, wide_mul_u16, wide_mul_u32};
pub mod utils;

/// Precomputed set of constants for each modulus which are used within
/// Montgomery multiplication + conversion.
///
/// Note: `n_prime` = (n^{-1}) \mod R, such that N N' = 1 \mod R
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
            n_prime: compute_inverse(modulus as u32, r),
        }
    }
}

/// Brings a value 0 \leq a < N into its "Montgomery form", aR \mod N.
/// We assume here that R = 2^{32} and that N < R.
///
/// The trick is to simply perform one wide mul by R^2 \mod N and then perform
/// a single reduction step, since now we are in the residue class a R^2 \mod N.
fn convert_to_montgomery_form_u16(a: u16, montgomery_parameters: &MontgomeryParametersU16) -> u16 {
    let a_times_r_squared = wide_mul_u16(a, montgomery_parameters.r_squared_mod_n);
    montgomery_reduction_u16(a_times_r_squared, montgomery_parameters)
}

const MOD_U32_BIT_MASK: u32 = (1 << 16) - 1;

fn mod_r_u16(a: u32) -> u16 {
    (a & MOD_U32_BIT_MASK)
        .try_into()
        .expect("Error: this should always work")
}

fn div_r_u16(a: u32) -> u16 {
    (a >> 16)
        .try_into()
        .expect("Error: this should always work")
}

/// Computes TR^{-1} \mod N. Note that R = 2^{32} in this case.
///
/// ## Assumptions
/// * 0 \leq T < RN
fn montgomery_reduction_u16(t_mont: u32, montgomery_parameters: &MontgomeryParametersU16) -> u16 {
    // (T + ((T \mod R N') \mod R) N) / R
    let t_mod_r = mod_r_u16(t_mont);
    // This implicitly happens mod R since overflow is silently truncated
    // Note: \alpha = ((T \mod R N') \mod R). We compute (T + \alpha N) / R
    let alpha = (t_mod_r as u32) * montgomery_parameters.n_prime;
    let numerator = (t_mont as u64) - wide_mul_u32(alpha, montgomery_parameters.modulus as u32);
    div_r_u16(numerator.try_into().expect("Yikes"))
}

fn main() {
    let modulus_digits = Generator::safe_prime(16).to_u32_digits();
    debug_assert_eq!(modulus_digits.len(), 1);
    let modulus = modulus_digits[0] as u16;

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

    dbg!(a);
    dbg!(b);
    dbg!(modulus);
    dbg!(ab_mod_modulus);
    dbg!(ab_mod_modulus_via_montgomery);

    assert_eq!(ab_mod_modulus, ab_mod_modulus_via_montgomery);
}
