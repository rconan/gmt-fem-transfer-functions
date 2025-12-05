//! Frequency response data products

#[cfg(feature = "faer")]
use faer::Mat;
#[cfg(feature = "nalgebra")]
use nalgebra::{ComplexField, DMatrix};
use num_complex::Complex;
use serde::Serialize;
use std::io::BufWriter;
use std::time::Instant;
use std::{env, f64, fmt::Display, fs::File, io, ops::Deref, path::Path};

use crate::{cli::Cli, structural::Structural};

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
    #[error("failed to write to Matlab data file")]
    Matlab(#[from] matio_rs::MatioError),
}

type Result<T> = std::result::Result<T, TransferFunctionDataError>;

/// Matrix and scale size interface
pub trait Dims {
    type D: std::fmt::Debug + Serialize;
    fn size(&self) -> Self::D;
}

#[cfg(feature = "faer")]
impl Dims for Mat<f64> {
    type D = (usize, usize);

    fn size(&self) -> Self::D {
        self.shape()
    }
}
#[cfg(feature = "nalgebra")]
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

#[cfg(feature = "faer")]
impl Cartesian2Polar for Mat<Complex<f64>> {
    type Output = Mat<f64>;

    fn magnitude(&self) -> Self::Output {
        let mut col_wise_data = self
            .col_iter()
            .flat_map(|col| col.iter().cloned().map(|c| c.norm()).collect::<Vec<_>>());
        let (nrows, ncols) = self.shape();
        Mat::from_fn(nrows, ncols, |_, _| col_wise_data.next().unwrap())
    }

    fn phase(&self) -> Self::Output {
        let mut col_wise_data = self
            .col_iter()
            .flat_map(|col| col.iter().map(|c| c.arg()).collect::<Vec<_>>());
        let (nrows, ncols) = self.shape();
        Mat::from_fn(nrows, ncols, |_, _| col_wise_data.next().unwrap())
    }
}

#[cfg(feature = "nalgebra")]
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
        self.norm()
    }

    fn phase(&self) -> Self::Output {
        self.arg()
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
#[derive(Debug, Serialize)]
pub struct FrequencyResponseVec<T: Cartesian2Polar>(
    #[serde(rename = "data")] Vec<FrequencyResponseData<T>>,
);
impl<T: Cartesian2Polar> Default for FrequencyResponseVec<T> {
    fn default() -> Self {
        Self(vec![])
    }
}

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
    fem_eigen_frequency_range: (f64, f64),
    #[cfg(feature = "nalgebra")]
    frequency_response: FrequencyResponseVec<DMatrix<Complex<f64>>>,
    #[cfg(feature = "faer")]
    frequency_response: FrequencyResponseVec<Mat<Complex<f64>>>,
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
        let now = Instant::now();
        match path.as_ref().extension() {
            Some(ext) if ext == "pkl" => {
                let file = File::create(&path)?;
                let mut buffer = BufWriter::new(file);
                serde_pickle::to_writer(&mut buffer, &self, Default::default())?;
            }
            #[cfg(feature = "nalgebra")]
            Some(ext) if ext == "mat" => self.dump_to_mat(&path)?,
            #[cfg(feature = "faer")]
            Some(ext) if ext == "mat" => unimplemented!(),
            Some(ext) => {
                return Err(TransferFunctionDataError::DataFileExtension(
                    ext.to_string_lossy().into_owned(),
                ));
            }
            None => return Err(TransferFunctionDataError::MissingFileExtension),
        };
        println!(
            "Frequency response written to {} in {}ms",
            path.as_ref().display(),
            now.elapsed().as_millis()
        );
        Ok(())
    }

    #[cfg(feature = "nalgebra")]
    pub fn dump_to_mat(self, path: impl AsRef<Path>) -> Result<()> {
        use matio_rs::{Mat, MatFile, MayBeFrom};
        let mut fields = vec![
            Mat::maybe_from("fem", self.fem)?,
            Mat::maybe_from("inputs", self.inputs)?,
            Mat::maybe_from("outputs", self.outputs)?,
            Mat::maybe_from("modal_damping_coefficient", self.modal_damping_coefficient)?,
            Mat::maybe_from("fem_eigen_frequency_range", self.fem_eigen_frequency_range)?,
        ];
        let mut data = vec![];
        for r in self.frequency_response.iter() {
            let data_fields = vec![
                Mat::maybe_from("frequency", r.frequency)?,
                Mat::maybe_from("magnitude", r.magnitude.clone())?,
                Mat::maybe_from("phase", r.phase.clone())?,
            ];
            data.push(Mat::maybe_from("data", data_fields)?);
        }
        let data_iter = Box::new(data.into_iter()) as Box<dyn Iterator<Item = Mat>>;
        fields.push(Mat::maybe_from("frequency_response", vec![data_iter])?);
        let mstruct = Mat::maybe_from("transfer_functions", fields)?;
        MatFile::save(path)?.write(mstruct);
        Ok(())
    }

    /// Adds the [frequency response](FrequencyResponseVec) to the data
    #[cfg(feature = "nalgebra")]
    pub fn add_response(
        self,
        frequency_response: FrequencyResponseVec<DMatrix<Complex<f64>>>,
    ) -> Self {
        Self {
            frequency_response,
            ..self
        }
    }
    #[cfg(feature = "faer")]
    pub fn add_response(self, frequency_response: FrequencyResponseVec<Mat<Complex<f64>>>) -> Self {
        Self {
            frequency_response,
            ..self
        }
    }

    /// Adds additional data from the structural model
    pub fn add_structural(self, structural: &Structural) -> Self {
        let c = 0.5 * f64::consts::FRAC_1_PI;
        Self {
            fem_eigen_frequency_range: (structural.w[0] * c, *structural.w.last().unwrap() * c),
            ..self
        }
    }
}
