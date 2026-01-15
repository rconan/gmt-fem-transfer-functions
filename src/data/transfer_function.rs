use std::{
    env, f64,
    fmt::Debug,
    fs::File,
    io::{self, BufWriter},
    path::Path,
    time::Instant,
};

use serde::Serialize;

use crate::{
    Cli, if64,
    structural::{Structural, StructuralError},
};

use super::{Cartesian2Polar, FrequencyResponseVec};

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

#[cfg(feature = "faer")]
type Matrix<T> = faer::Mat<T>;
#[cfg(feature = "nalgebra")]
type Matrix<T> = nalgebra::DMatrix<T>;

#[derive(Debug, Serialize)]

pub struct ModalMatrix {
    pub(crate) mat: Matrix<f64>,
    pub(crate) nodes: Vec<f64>,
}
impl ModalMatrix {
    pub fn new(mat: Matrix<f64>, nodes: Vec<f64>) -> Self {
        Self { mat, nodes }
    }
}

/// GMT FEM transfer function data export
#[derive(Debug, Default, Serialize)]
pub struct TransferFunctionData<T>
where
    T: Clone + PartialEq + Debug + 'static,
    Matrix<T>: Cartesian2Polar,
{
    fem: String,
    inputs: Vec<String>,
    outputs: Vec<String>,
    modal_damping_coefficient: f64,
    fem_eigen_frequency_range: (f64, f64),
    frequency_response: FrequencyResponseVec<Matrix<T>>,
    pub(crate) b: Option<ModalMatrix>,
    // modal displacements matrix
    pub(crate) c: Option<ModalMatrix>,
}

impl<T> From<&Cli> for TransferFunctionData<T>
where
    T: Clone + PartialEq + Default + Debug,
    Matrix<T>: Cartesian2Polar,
{
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

impl TransferFunctionData<if64> {
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
            Some(ext) if ext == "mat" => self.dump_to_mat(&path)?,
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
            let data_fields = if let Some(phase) = &r.phase {
                vec![
                    Mat::maybe_from("frequency", r.frequency)?,
                    Mat::maybe_from("magnitude", &r.magnitude)?,
                    Mat::maybe_from("phase", phase)?,
                ]
            } else {
                vec![
                    Mat::maybe_from("frequency", r.frequency)?,
                    Mat::maybe_from("magnitude", &r.magnitude)?,
                ]
            };
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
    pub fn add_response(self, frequency_response: FrequencyResponseVec<nalgebra::DMatrix<if64>>) -> Self {
        Self {
            frequency_response,
            ..self
        }
    }
    #[cfg(feature = "faer")]
    pub fn add_response(self, frequency_response: FrequencyResponseVec<faer::Mat<if64>>) -> Self {
        Self {
            frequency_response,
            ..self
        }
    }

    /// Adds additional data from the structural model
    pub fn add_structural(
        self,
        structural: &Structural,
        b: bool,
        c: bool,
    ) -> std::result::Result<Self, StructuralError> {
        let sc = 0.5 * f64::consts::FRAC_1_PI;
        Ok(Self {
            fem_eigen_frequency_range: (structural.w[0] * sc, *structural.w.last().unwrap() * sc),
            b: b.then(|| {
                Ok::<_, StructuralError>(ModalMatrix::new(
                    structural.force_to_mode(),
                    structural.inputs_nodes()?,
                ))
            })
            .transpose()?,
            c: c.then(|| {
                Ok::<_, StructuralError>(ModalMatrix::new(
                    structural.mode_to_displacement(),
                    structural.outputs_nodes()?,
                ))
            })
            .transpose()?,
            ..self
        })
    }
}

#[cfg(feature = "faer")]
impl TransferFunctionData<f64> {
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
            Some(ext) if ext == "mat" => self.dump_to_mat(&path)?,
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
            let data_fields = if let Some(phase) = &r.phase {
                vec![
                    Mat::maybe_from("frequency", r.frequency)?,
                    Mat::maybe_from("magnitude", &r.magnitude)?,
                    Mat::maybe_from("phase", phase)?,
                ]
            } else {
                vec![
                    Mat::maybe_from("frequency", r.frequency)?,
                    Mat::maybe_from("magnitude", &r.magnitude)?,
                ]
            };
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
    pub fn add_response(self, frequency_response: FrequencyResponseVec<nalgebra::DMatrix<f64>>) -> Self {
        Self {
            frequency_response,
            ..self
        }
    }
    #[cfg(feature = "faer")]
    pub fn add_response(self, frequency_response: FrequencyResponseVec<faer::Mat<f64>>) -> Self {
        Self {
            frequency_response,
            ..self
        }
    }

    /// Adds additional data from the structural model
    pub fn add_structural(
        self,
        structural: &Structural,
        b: bool,
        c: bool,
    ) -> std::result::Result<Self, StructuralError> {
        let sc = 0.5 * f64::consts::FRAC_1_PI;
        Ok(Self {
            fem_eigen_frequency_range: (structural.w[0] * sc, *structural.w.last().unwrap() * sc),
            b: b.then(|| {
                Ok::<_, StructuralError>(ModalMatrix::new(
                    structural.force_to_mode(),
                    structural.inputs_nodes()?,
                ))
            })
            .transpose()?,
            c: c.then(|| {
                Ok::<_, StructuralError>(ModalMatrix::new(
                    structural.mode_to_displacement(),
                    structural.outputs_nodes()?,
                ))
            })
            .transpose()?,
            ..self
        })
    }
}
