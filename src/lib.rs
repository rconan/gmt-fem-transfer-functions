//! GMT FEM frequency response

use clap::ValueEnum;

pub mod cli;
#[doc(hidden)]
pub use cli::Cli;
pub mod data;
pub mod frequency_response;
pub mod structural;

include!(concat!(env!("OUT_DIR"), "/fem_io.rs"));

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("failed to process command line interface arguments")]
    Cli(#[from] cli::CliError),
    #[error("failed to build structural model")]
    Structural(#[from] structural::StructuralError),
}
