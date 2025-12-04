//! Command line interface

use std::{io, time::Instant};

use crate::{Inputs, Outputs, frequency_response::Frequencies};
use clap::Parser;
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

#[derive(Debug, thiserror::Error)]
pub enum CliError {
    #[error("failed to load optical sensitivities")]
    LoadOpticalSensitivities(#[from] io::Error),
    #[error("failed to decompress optical sensitivities")]
    DecompressOpticalSensitivities(#[from] lz4_flex::block::DecompressError),
    #[error("failed to deserialize optical sensitivities")]
    DeserOpticalSensitivities(#[from] bincode::error::DecodeError),
    #[error("")]
    Lom,
}

#[derive(Default, Debug, Serialize, Deserialize)]
pub struct Lom {
    pub tip_tilt: Vec<f64>,
    pub segment_tip_tilt: Vec<f64>,
    pub segment_piston: Vec<f64>,
}
impl Lom {
    pub fn new() -> Result<Self, CliError> {
        let now = Instant::now();
        let buffer: &[u8] = include_bytes!("lom.lz4");
        let decompressed = lz4_flex::decompress_size_prepended(buffer)?;
        let (lom, _): (Lom, usize) =
            bincode::serde::decode_from_slice(&decompressed, bincode::config::standard())?;
        println!(
            "loaded linear optical sensitivity matrices in {}mus",
            now.elapsed().as_micros()
        );
        Ok(lom)
    }
}

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
    pub fn lom_sensitivies(&self) -> Result<Option<DMatrix<f64>>, CliError> {
        let mut lom = Option::<Lom>::None;
        let mats: Vec<DMatrix<f64>> = self
            .outputs
            .iter()
            .map(|output| match output {
                Outputs::TipTilt => Ok(Some(DMatrix::<f64>::from_column_slice(
                    2,
                    84,
                    &{
                        if lom.is_none() {
                            lom = Some(Lom::new()?);
                            lom.as_ref()
                        } else {
                            lom.as_ref()
                        }
                    }
                    .unwrap()
                    .tip_tilt,
                ))),
                Outputs::SegmentTipTilt => Ok(Some(DMatrix::<f64>::from_column_slice(
                    14,
                    84,
                    &{
                        if lom.is_none() {
                            lom = Some(Lom::new()?);
                            lom.as_ref()
                        } else {
                            lom.as_ref()
                        }
                    }
                    .unwrap()
                    .segment_tip_tilt,
                ))),
                Outputs::SegmentPiston => Ok(Some(DMatrix::<f64>::from_column_slice(
                    7,
                    84,
                    &{
                        if lom.is_none() {
                            lom = Some(Lom::new()?);
                            lom.as_ref()
                        } else {
                            lom.as_ref()
                        }
                    }
                    .unwrap()
                    .segment_piston,
                ))),
                _ => Ok(None),
            })
            .filter_map(|mat| mat.transpose())
            .collect::<Result<Vec<_>, CliError>>()?;
        Ok(lom.map(|_| {
            let rows: Vec<_> = mats
                .into_iter()
                .flat_map(|mat| mat.row_iter().map(|r| r.into_owned()).collect::<Vec<_>>())
                .collect();
            DMatrix::<f64>::from_rows(&rows)
        }))
    }
}
