//! GMT FEM frequency response

use clap::ValueEnum;

mod cli;
pub use cli::Cli;
pub mod data;
pub mod frequency_response;
pub mod structural;

include!(concat!(env!("OUT_DIR"), "/fem_io.rs"));
