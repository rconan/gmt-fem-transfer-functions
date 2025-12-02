use crate::{Inputs, Outputs, frequency_response::Frequencies};
use clap::Parser;

#[derive(Parser)]
pub struct Cli {
    /// FEM inputs
    #[arg(short, long)]
    pub inputs: Vec<Inputs>,
    /// FEM outputs
    #[arg(short, long)]
    pub outputs: Vec<Outputs>,
    /// FEM modal damping coeffcient
    #[arg(short = 'z', long, default_value_t = 0.02f64)]
    pub structural_damping: f64,
    /// Frequencies [Hz]
    #[command(subcommand)]
    pub frequencies: Frequencies,
    /// data file, either a Matlab (.mat) or Python pickle (.pkl) file
    #[arg(short, long, default_value_t = String::from("gmt_frequency_response.pkl"))]
    pub filename: String,
}
