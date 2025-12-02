use nalgebra::{Complex, ComplexField, DMatrix};
use serde::Serialize;
use std::{env, fmt::Display, fs::File, io, path::Path};

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

pub trait Cartesian2Polar {
    type Output: std::fmt::Debug + Serialize;
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

#[derive(Debug, Serialize)]
pub struct FrequencyResponseData<T: Cartesian2Polar> {
    frequency: f64,
    magnitude: <T as Cartesian2Polar>::Output,
    phase: <T as Cartesian2Polar>::Output,
}
impl<T: Cartesian2Polar> FrequencyResponseData<T> {
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

#[derive(Debug, Default, Serialize)]
pub struct TransferFunctionData {
    fem: String,
    inputs: Vec<String>,
    outputs: Vec<String>,
    modal_damping_coefficient: f64,
    frequency_response: Vec<FrequencyResponseData<DMatrix<Complex<f64>>>>,
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

    pub fn add_response(
        self,
        frequency_response: Vec<FrequencyResponseData<DMatrix<Complex<f64>>>>,
    ) -> Self {
        Self {
            frequency_response,
            ..self
        }
    }
}
