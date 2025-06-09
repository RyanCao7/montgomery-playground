//! The idea behind Barrett is to compute ab \mod N, but to do so without
//! needing to divide by N.
//!
//! We instead compute \lfloor R / N \rfloor, where R = 2^{32} or whatever our
//! modulus range is supposed to be in.

/// Implicitly, our "radix" (`c` above) is 2^{16}.
pub struct BarrettParamsU16 {
    modulus: u16,
    mu: u32,
}
const LHS_SHIFT_U16: u32 = 15;
const RHS_SHIFT_U16: u32 = 17;

pub fn barrett_reduction_u16(val: u32, barrett_params: &BarrettParamsU16) -> u16 {
    let lhs_term: u16 = (val >> LHS_SHIFT_U16).try_into().expect("This should work");
    let rhs_term: u16 = (barrett_params.mu >> RHS_SHIFT_U16)
        .try_into()
        .expect("This should work too");
    todo!()
}
