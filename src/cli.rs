//! Command line interface

use crate::{Inputs, Outputs, frequency_response::Frequencies};
use clap::Parser;

/// Command line interface
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
    /// FEM minimum eigen frequency (lower frequencies are dropped)
    #[arg(long)]
    pub eigen_frequency_min: Option<f64>,
    /// FEM maximum eigen frequency (higher frequencies are truncated)
    #[arg(long)]
    pub eigen_frequency_max: Option<f64>,
    /// Frequencies \[Hz\]
    #[command(subcommand)]
    pub frequencies: Frequencies,
    /// data file, either a Matlab (.mat) or Python pickle (.pkl) file
    #[arg(short, long, default_value_t = String::from("gmt_frequency_response.pkl"))]
    pub filename: String,
}

impl Cli {
    /// Returns the names of the FEM inputs
    pub fn fem_inputs(&self) -> Vec<String> {
        self.inputs.iter().map(|io| io.name()).collect()
    }
    /// Returns the names of the FEM outputs
    pub fn fem_outputs(&self) -> Vec<String> {
        let mut outs: Vec<_> = self
            .outputs
            .iter()
            .map(|io| io.name())
            .filter(|name| {
                !(name == "tip-tilt" || name == "segment_tip-tilt" || name == "segment_piston")
            })
            .collect();
        #[cfg(fem)]
        for output in &self.outputs {
            match output {
                Outputs::TipTilt | Outputs::SegmentTipTilt | Outputs::SegmentPiston => {
                    let var = Outputs::OSSM1Lcl;
                    outs.push(var.name());
                    let var = Outputs::MCM2Lcl6D;
                    outs.push(var.name());
                }
                _ => {}
            }
        }
        outs.dedup();
        outs
    }
    /// Returns the names of the linear optical model outputs
    pub fn lom_outputs(&self) -> Vec<String> {
        self.outputs
            .iter()
            .map(|io| io.name())
            .filter(|name| {
                name == "tip-tilt" || name == "segment_tip-tilt" || name == "segment_piston"
            })
            .collect()
    }
}
