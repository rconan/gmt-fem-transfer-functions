//! Frequency response functionalities

use indicatif::{ParallelProgressIterator, ProgressStyle};
use rayon::prelude::*;
use std::{f64::consts::PI, ops::Mul};

use crate::{
    data::{Cartesian2Polar, FrequencyResponseData, FrequencyResponseVec, Get},
    if64,
};

mod frequencies;
pub use frequencies::Frequencies;

const DPI: f64 = 2f64 * PI;

pub trait JOmega {
    /// Transfer function type
    type Output;

    /// Returns the frequency response
    ///
    /// The argument is the imaginary frequency in radians
    fn j_omega(&self, jw: if64) -> Self::Output;
}

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
        let style = ProgressStyle::with_template("|{bar} {pos}|")
            .unwrap()
            .progress_chars("-.-");
        let data = match frequencies {
            Frequencies::Single { value: nu } => {
                let jw = if64::new(0f64, DPI * nu);
                vec![FrequencyResponseData::new(nu, self.j_omega(jw))]
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
                        FrequencyResponseData::new(nu, self.j_omega(jw))
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
                        FrequencyResponseData::new(nu, self.j_omega(jw))
                    })
                    .collect()
            }
            Frequencies::Set { values: nu } => nu
                .into_par_iter()
                .progress_with_style(style)
                .map(|nu| {
                    let jw = if64::new(0f64, DPI * nu);
                    FrequencyResponseData::new(nu, self.j_omega(jw))
                })
                .collect(),
            _ => panic!("frequencies not set, aborting"),
        };
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
impl<T: JOmegaSvd> FrequencyResponseSvd for T {}

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

#[cfg(test)]
mod tests {
    // use std::fs::File;

    use super::*;

    #[test]
    fn folp_tf() {
        let folp = FirstOrderLowPass::new();

        let tf = folp.frequency_response(Frequencies::logspace(1., 8e3, 1000));

        // let mut file = File::create("folp_tf.pkl").unwrap();
        // serde_pickle::to_writer(&mut file, &(nu, tf), Default::default()).unwrap();
    }

    #[test]
    fn bessel_tf() {
        let bessel = BesselFilter::new();

        let tf = bessel.frequency_response(Frequencies::logspace(1., 8e3, 1000));

        // let mut file = File::create("bessel_tf.pkl").unwrap();
        // serde_pickle::to_writer(&mut file, &(nu, tf), Default::default()).unwrap();
    }

    #[test]
    fn pic_tf() {
        let pic = PICompensator::new();

        let tf = pic.frequency_response(Frequencies::logspace(1., 8e3, 1000));

        // let mut file = File::create("pic_tf.pkl").unwrap();
        // serde_pickle::to_writer(&mut file, &(nu, tf), Default::default()).unwrap();
    }
}
