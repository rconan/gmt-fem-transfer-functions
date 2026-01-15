//! Frequency response functionalities

use std::{f64::consts::PI, ops::Mul};

use crate::{
    data::{Cartesian2Polar, FrequencyResponseData, FrequencyResponseVec, Get},
    if64,
};

mod frequencies;
pub use frequencies::Frequencies;
mod svd;
pub use svd::{FrequencyResponseSvd, JOmegaSvd};

const DPI: f64 = 2f64 * PI;

pub trait JOmega {
    /// Transfer function type
    type Output;

    /// Returns the frequency response
    ///
    /// The argument is the imaginary frequency in radians
    fn j_omega(&self, jw: if64) -> Self::Output;
}

/// Frequency response interface definition
pub trait FrequencyResponse: JOmega {
    /// Returns the frequencies and the frequency response
    ///
    /// The argument is frequencies in Hz
    fn frequency_response<T: Into<Frequencies>>(&self, nu: T) -> FrequencyResponseVec<Self::Output>
    where
        <Self as JOmega>::Output: Cartesian2Polar + Send,
        <<Self as JOmega>::Output as Cartesian2Polar>::Output: Get + Send,
        Self: Sync,
    {
        let frequencies: Frequencies = nu.into();
        let data = frequencies.data(|nu, jw| FrequencyResponseData::new(nu, self.j_omega(jw)));
        FrequencyResponseVec::new(data)
    }

    /// Returns the first derivation of the frequency response
    fn j_omega_first(&self, jw: if64) -> <<Self as JOmega>::Output as Mul<if64>>::Output
    where
        <Self as JOmega>::Output: Mul<if64>,
    {
        self.j_omega(jw) * jw
    }
    /// Returns the second derivation of the frequency response
    fn j_omega_second(
        &self,
        jw: if64,
    ) -> <<<Self as JOmega>::Output as Mul<if64>>::Output as Mul<if64>>::Output
    where
        <Self as JOmega>::Output: Mul<if64>,
        <<Self as JOmega>::Output as Mul<if64>>::Output: Mul<if64>,
    {
        self.j_omega_first(jw) * jw
    }
}
impl<T: JOmega> FrequencyResponse for T {}
