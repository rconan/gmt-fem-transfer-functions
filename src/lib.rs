//! GMT FEM frequency response

use clap::ValueEnum;

pub mod cli;
#[doc(inline)]
pub use cli::Cli;
pub mod data;
pub mod frequency_response;
pub mod structural;

include!(concat!(env!("OUT_DIR"), "/fem_io.rs"));


#[allow(non_camel_case_types)]
pub type if64 = num_complex::Complex<f64>;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("failed to process command line interface arguments")]
    Cli(#[from] cli::CliError),
    #[error("failed to build structural model")]
    Structural(#[from] structural::StructuralError),
}

impl TryFrom<&Cli> for structural::Structural {
    type Error = crate::Error;

    fn try_from(args: &Cli) -> Result<Self, Self::Error> {
        Ok(structural::StructuralBuilder {
            built: structural::Structural {
                inputs: args.fem_inputs(),
                outputs: args.fem_outputs(),
                z: args.structural_damping,
                optical_senses: args.lom_sensitivies()?,
                ..Default::default()
            },
            min_eigen_frequency: args.eigen_frequency_min,
            max_eigen_frequency: args.eigen_frequency_max,
            ..Default::default()
        }
        .build()?)
    }
}

