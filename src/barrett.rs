//! The idea behind Barrett is to compute ab \mod N, but to do so without
//! needing to divide by N.
//!
//! We instead compute \lfloor R / N \rfloor, where R = 2^{32} or whatever our
//! modulus range is supposed to be in.

use crate::utils::{div_r_u32, wide_mul_u32};
use std::u32::MAX as U32_MAX;

pub struct BarrettParamsU16 {
    modulus: u16,
    mu: u32, // This is \floor{R / N}
}
impl BarrettParamsU16 {
    pub fn new(modulus: u16) -> Self {
        Self {
            modulus,
            mu: U32_MAX / (modulus as u32), // This is always okay since `modulus` is prime
        }
    }
}

/// Basic reduction formula (we are not implementing the "advanced" one because
/// I don't see how it can work with the most pessimistic ranges):
///
/// x - \floor[ (x * \mu) / R ] * N
pub fn barrett_reduction_u16(val: u32, barrett_params: &BarrettParamsU16) -> u16 {
    let numerator = wide_mul_u32(val, barrett_params.mu);
    let numerator_divided_by_r = div_r_u32(numerator);
    let product_with_n = numerator_divided_by_r * (barrett_params.modulus as u32);
    (val - product_with_n) as u16 // No conditional subtraction needed since R > N^2
}
