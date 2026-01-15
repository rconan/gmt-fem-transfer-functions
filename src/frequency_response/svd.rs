use crate::{
    data::{Cartesian2Polar, FrequencyResponseData, FrequencyResponseVec, Get},
    if64,
};

use super::{Frequencies, JOmega};

pub trait JOmegaSvd: JOmega {
    /// SVD matrices type
    type Svd;

    /// Returns the frequency response singular values
    ///
    /// The argument is the imaginary frequency in radians
    fn j_omega_svd(
        &self,
        jw: if64,
        u: bool,
        v: bool,
    ) -> (Self::Svd, Option<Self::Svd>, Option<Self::Svd>);
}

pub trait FrequencyResponseSvd: JOmegaSvd {
    fn frequency_response_svd<T: Into<Frequencies>>(
        &self,
        nu: T,
        u: bool,
        v: bool,
    ) -> FrequencyResponseVec<Self::Svd>
    where
        <Self as JOmegaSvd>::Svd: Cartesian2Polar + Send,
        <<Self as JOmegaSvd>::Svd as Cartesian2Polar>::Output: Get + Send,
        Self: Sync,
    {
        let frequencies: Frequencies = nu.into();
        let data = frequencies
            .data(|nu, jw| FrequencyResponseData::new_svd(nu, self.j_omega_svd(jw, u, v)));
        FrequencyResponseVec::new(data)
    }
}

impl<T: JOmegaSvd> FrequencyResponseSvd for T {}
