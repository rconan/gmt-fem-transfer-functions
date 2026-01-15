use std::f64::consts::PI;

use indicatif::{ParallelProgressIterator, ProgressStyle};
use rayon::prelude::*;

use crate::{
    data::{Cartesian2Polar, FrequencyResponseData, FrequencyResponseVec, Get},
    if64,
};

pub use super::{Frequencies, JOmega};

const DPI: f64 = 2f64 * PI;

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
        let style = ProgressStyle::with_template("|{bar} {pos}|")
            .unwrap()
            .progress_chars("-.-");
        let data = match frequencies {
            Frequencies::Single { value: nu } => {
                let jw = if64::new(0f64, DPI * nu);
                vec![FrequencyResponseData::new_svd(
                    nu,
                    self.j_omega_svd(jw, u, v),
                )]
            }
            Frequencies::LogSpace { lower, upper, n } => {
                assert!(upper > lower);
                let log_step = (upper.log10() - lower.log10()) / (n - 1) as f64;
                (0..n)
                    .into_par_iter()
                    .progress_with_style(style)
                    .map(|i| {
                        let log_nu = lower.log10() + log_step * i as f64;
                        let nu = 10f64.powf(log_nu);
                        let jw = if64::new(0f64, DPI * nu);
                        FrequencyResponseData::new_svd(nu, self.j_omega_svd(jw, u, v))
                    })
                    .collect()
            }
            Frequencies::LinSpace { lower, upper, n } => {
                assert!(upper > lower);
                let step = (upper - lower) / (n - 1) as f64;
                (0..n)
                    .into_par_iter()
                    .progress_with_style(style)
                    .map(|i| {
                        let nu = lower + step * i as f64;
                        let jw = if64::new(0f64, DPI * nu);
                        FrequencyResponseData::new_svd(nu, self.j_omega_svd(jw, u, v))
                    })
                    .collect()
            }
            Frequencies::Set { values: nu } => nu
                .into_par_iter()
                .progress_with_style(style)
                .map(|nu| {
                    let jw = if64::new(0f64, DPI * nu);
                    FrequencyResponseData::new_svd(nu, self.j_omega_svd(jw, u, v))
                })
                .collect(),
            _ => panic!("frequencies not set, aborting"),
        };
        FrequencyResponseVec::new(data)
    }
}

impl<T: JOmegaSvd> FrequencyResponseSvd for T {}
