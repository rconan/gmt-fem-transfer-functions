//! Command line interface

use crate::{
    Inputs, Outputs,
    frequency_response::Frequencies,
    structural::{Structural, StructuralBuilder, StructuralError},
};
use clap::Parser;
use gmt_lom::{
    LinearOpticalModelError, Loader, LoaderTrait, OpticalSensitivities, OpticalSensitivity,
};
use nalgebra::{DMatrix, RowDVector};

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
            if let Outputs::TipTilt | Outputs::SegmentTipTilt | Outputs::SegmentPiston = output {
                let var = Outputs::OSSM1Lcl;
                outs.push(var.name());
                let var = Outputs::MCM2Lcl6D;
                outs.push(var.name());
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
    pub fn lom_sensitivies(&self) -> Result<Option<DMatrix<f64>>, LinearOpticalModelError> {
        let sensitivities = Loader::<OpticalSensitivities>::default().load()?;
        let rows: Vec<_> = self
            .outputs
            .iter()
            .filter_map(|output| match output {
                Outputs::TipTilt => {
                    Some(sensitivities[OpticalSensitivity::TipTilt(vec![])].clone())
                }
                Outputs::SegmentTipTilt => {
                    Some(sensitivities[OpticalSensitivity::SegmentTipTilt(vec![])].clone())
                }
                Outputs::SegmentPiston => {
                    Some(sensitivities[OpticalSensitivity::SegmentPiston(vec![])].clone())
                }
                _ => None,
            })
            .flat_map(|sens| match sens {
                // gmt_lom::OpticalSensitivity::Wavefront(items) => items.,
                OpticalSensitivity::TipTilt(items) => items
                    .chunks(84)
                    .map(|x| RowDVector::from_row_slice(x))
                    .collect::<Vec<_>>(),
                OpticalSensitivity::SegmentTipTilt(items) => items
                    .chunks(84)
                    .map(|x| RowDVector::from_row_slice(x))
                    .collect::<Vec<_>>(),
                OpticalSensitivity::SegmentPiston(items) => items
                    .chunks(84)
                    .map(|x| RowDVector::from_row_slice(x))
                    .collect::<Vec<_>>(),
                _ => vec![],
            })
            .collect();
        Ok((!rows.is_empty()).then(|| DMatrix::<f64>::from_rows(rows.as_slice())))
    }
}

impl TryFrom<&Cli> for Structural {
    type Error = StructuralError;

    fn try_from(args: &Cli) -> Result<Self, Self::Error> {
        Ok(StructuralBuilder {
            built: Structural {
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
