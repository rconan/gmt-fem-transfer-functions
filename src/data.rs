//! Frequency response data products

use nalgebra::{Complex, ComplexField, DMatrix};
use serde::Serialize;
use std::{env, fmt::Display, fs::File, io, ops::Deref, path::Path};

use crate::cli::Cli;

#[derive(Debug, thiserror::Error)]
pub enum TransferFunctionDataError {
    #[error(r#"found data file extension: "{0}", expected "mat" or "pkl""#)]
    DataFileExtension(String),
    #[error(r#"missing data file extension: "mat" or "pkl""#)]
    MissingFileExtension,
    #[error("failed to create data file: {0}")]
    CreateDataFile(#[from] io::Error),
    #[error("failed to serialize data to pickle file")]
    SerPkl(#[from] serde_pickle::Error),
}

type Result<T> = std::result::Result<T, TransferFunctionDataError>;

/// Matrix and scale size interface
pub trait Dims {
    type D: std::fmt::Debug + Serialize;
    fn size(&self) -> Self::D;
}

impl Dims for DMatrix<f64> {
    type D = (usize, usize);

    fn size(&self) -> Self::D {
        self.shape()
    }
}

impl Dims for f64 {
    type D = usize;

    fn size(&self) -> Self::D {
        1
    }
}

/// Cartesian to polar transformation interface
pub trait Cartesian2Polar {
    type Output: Dims + std::fmt::Debug + Serialize;
    fn magnitude(&self) -> Self::Output;
    fn phase(&self) -> Self::Output;
}

impl Cartesian2Polar for DMatrix<Complex<f64>> {
    type Output = DMatrix<f64>;

    fn magnitude(&self) -> Self::Output {
        self.map(|x| x.modulus())
    }

    fn phase(&self) -> Self::Output {
        self.map(|x| x.argument())
    }
}

impl Cartesian2Polar for Complex<f64> {
    type Output = f64;

    fn magnitude(&self) -> Self::Output {
        self.modulus()
    }

    fn phase(&self) -> Self::Output {
        self.argument()
    }
}

/// Frequency response data point
///
/// Frequency response magnitude and phase matrices at one frequency
#[derive(Debug, Serialize)]
pub struct FrequencyResponseData<T: Cartesian2Polar> {
    frequency: f64,
    magnitude: <T as Cartesian2Polar>::Output,
    phase: <T as Cartesian2Polar>::Output,
}
impl<T: Cartesian2Polar> FrequencyResponseData<T> {
    /// Creates a [FrequencyResponseData] instance from a frequency and response complex matrix
    pub fn new(frequency: f64, response: T) -> Self {
        Self {
            frequency,
            magnitude: response.magnitude(),
            phase: response.phase(),
        }
    }
}
impl<T> Display for FrequencyResponseData<T>
where
    T: Cartesian2Polar,
    <T as Cartesian2Polar>::Output: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{},{},{}", self.frequency, self.magnitude, self.phase)
    }
}

/// Collection of [FrequencyResponseData]
#[derive(Debug, Default, Serialize)]
pub struct FrequencyResponseVec<T: Cartesian2Polar>(Vec<FrequencyResponseData<T>>);

impl<T: Cartesian2Polar> FrequencyResponseVec<T> {
    /// Creates a new [FrequencyResponseVec] instance from a vector of [FrequencyResponseData]
    pub fn new(frequency_response_datas: Vec<FrequencyResponseData<T>>) -> Self {
        Self(frequency_response_datas)
    }
    pub fn frequencies(&self) -> Vec<f64> {
        self.iter().map(|fr| fr.frequency).collect()
    }
}

impl<T: Cartesian2Polar> Deref for FrequencyResponseVec<T> {
    type Target = [FrequencyResponseData<T>];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Cartesian2Polar> FromIterator<FrequencyResponseData<T>> for FrequencyResponseVec<T> {
    fn from_iter<I: IntoIterator<Item = FrequencyResponseData<T>>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<T: Cartesian2Polar> Display for FrequencyResponseVec<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "GMT FEM frequency response matrix {}x{:?}",
            self.len(),
            self[0].magnitude.size()
        )?;
        match self.len() {
            n if n == 1 => {
                writeln!(f, " @ {:.2}Hz", self[0].frequency)
            }
            n if n < 6 => {
                writeln!(f, " @ {:.2?}Hz", self.frequencies())
            }
            _ => {
                writeln!(
                    f,
                    " @ [{:.2},{:.2}]Hz",
                    self[0].frequency,
                    self.last().unwrap().frequency
                )
            }
        }
    }
}

/// GMT FEM transfer function data export
#[derive(Debug, Default, Serialize)]
pub struct TransferFunctionData {
    fem: String,
    inputs: Vec<String>,
    outputs: Vec<String>,
    modal_damping_coefficient: f64,
    frequency_response: FrequencyResponseVec<DMatrix<Complex<f64>>>,
}

impl From<&Cli> for TransferFunctionData {
    fn from(args: &Cli) -> Self {
        let fem = {
            let fem_repo = env::var("FEM_REPO").unwrap();
            let fem_path = Path::new(&fem_repo);
            fem_path.file_name().unwrap().to_string_lossy().into_owned()
        };
        let inputs: Vec<_> = args.inputs.iter().map(|x| x.name()).collect();
        let outputs: Vec<_> = args.outputs.iter().map(|x| x.name()).collect();
        Self {
            fem,
            inputs,
            outputs,
            modal_damping_coefficient: args.structural_damping,
            ..Default::default()
        }
    }
}

impl TransferFunctionData {
    /// Writes the date to either a pickle or matlab file
    ///
    /// The file extension, "pkl" or "mat", sets the file type
    pub fn dump(self, path: impl AsRef<Path>) -> Result<()> {
        match path.as_ref().extension() {
            Some(ext) if ext == "pkl" => {
                let mut file = File::create(path)?;
                serde_pickle::to_writer(&mut file, &self, Default::default())?;
                Ok(())
            }
            Some(ext) if ext == "mat" => todo!(),
            Some(ext) => Err(TransferFunctionDataError::DataFileExtension(
                ext.to_string_lossy().into_owned(),
            )),
            None => Err(TransferFunctionDataError::MissingFileExtension),
        }
    }

    /// Adds the [frequency response](FrequencyResponseVec) to the data
    pub fn add_response(
        self,
        frequency_response: FrequencyResponseVec<DMatrix<Complex<f64>>>,
    ) -> Self {
        Self {
            frequency_response,
            ..self
        }
    }
}
