use crate::utils::{MOD_U32_BIT_MASK, compute_inverse, div_r_u16, mod_r_u16, wide_mul_u16};

/// Precomputed set of constants for each modulus which are used within
/// Montgomery multiplication + conversion.
///
/// Note: `n_prime` = -(n^{-1}) \mod R, such that N N' = -1 \mod R
pub struct MontgomeryParametersU16 {
    modulus: u16,
    r_squared_mod_n: u16,
    n_prime: u32,
}

impl MontgomeryParametersU16 {
    pub fn new(modulus: u16) -> Self {
        let r = 1_u32 << 16;
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
pub fn convert_to_montgomery_form_u16(
    a: u16,
    montgomery_parameters: &MontgomeryParametersU16,
) -> u16 {
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
pub fn montgomery_reduction_u16(
    t_mont: u32,
    montgomery_parameters: &MontgomeryParametersU16,
) -> u16 {
    debug_assert!(t_mont < (1_u32 << 16) * (montgomery_parameters.modulus as u32));
    let t_mod_r = mod_r_u16(t_mont);
    // Note: \alpha = ((T \mod R N') \mod R).
    let alpha = mod_r_u16((t_mod_r as u32) * montgomery_parameters.n_prime);
    // We compute (T + \alpha N) / R
    let numerator = t_mont + wide_mul_u16(alpha, montgomery_parameters.modulus);
    // The lower 16 bits should all be zero, i.e. divisible by 2^{16}
    debug_assert!(numerator & MOD_U32_BIT_MASK == 0);
    div_r_u16(numerator) // Finally, right-shift by R
}
