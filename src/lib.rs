use clap::ValueEnum;

pub mod cli;
pub mod data;
pub mod frequency_response;
pub mod structural;

pub trait BuilderTrait {
    /// Sets the FEM modal damping coefficient
    fn damping(self, z: f64) -> Self;
    /// Truncates the eigen frequencies to and including `max_eigen_frequency`
    ///
    /// The number of modes is set accordingly
    fn max_eigen_frequency(self, max_eigen_frequency: Option<f64>) -> Self;
    /// Drops the eigen frequencies less than `min_eigen_frequency`
    ///
    /// The number of modes is set accordingly
    fn min_eigen_frequency(self, min_eigen_frequency: Option<f64>) -> Self;
    /// Sets the filename where [Structural] is seralize to
    fn filename<S: Into<String>>(self, file_name: S) -> Self;
    /// Enables the compensation of the static gain mismatch
    ///
    /// An optional delay [s] may be added
    fn enable_static_gain_mismatch_compensation(self, maybe_delay: Option<f64>) -> Self;
}

include!(concat!(env!("OUT_DIR"), "/fem_io.rs"));
