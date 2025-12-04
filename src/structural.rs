//! FEM structural dynamic model

use std::{f64::consts, fmt::Display};

use gmt_dos_clients_fem::{Model, Switch};
use gmt_fem::FEM;
use nalgebra::{DMatrix, DMatrixView};
use num_complex::Complex;
use serde::{Deserialize, Serialize};

use crate::frequency_response::{FrequencyResponse, if64};

#[derive(Debug, thiserror::Error)]
pub enum StructuralError {
    #[error(transparent)]
    FEM(#[from] gmt_fem::FemError),
    // #[error(transparent)]
    // Decode(#[from] bincode::error::DecodeError),
    // #[error(transparent)]
    // Encode(#[from] bincode::error::EncodeError),
    #[error(transparent)]
    IO(#[from] std::io::Error),
    #[error("inputs and outputs do not match model in {0}")]
    IOMismatch(String),
}
type Result<T> = std::result::Result<T, StructuralError>;

#[derive(Debug, Deserialize, Serialize)]
pub struct StaticGainCompensation {
    pub(crate) delay: Option<f64>,
    pub(crate) delta_gain: DMatrix<if64>,
}
impl Default for StaticGainCompensation {
    fn default() -> Self {
        Self {
            delay: Default::default(),
            delta_gain: DMatrix::<if64>::zeros(1, 1),
        }
    }
}

/// FEM structural dynamic model
#[derive(Debug, Default, Deserialize, Serialize)]
pub struct Structural {
    // inputs labels
    pub(crate) inputs: Vec<String>,
    // outputs labels
    pub(crate) outputs: Vec<String>,
    // modal forces matrix
    pub(crate) b: DMatrix<if64>,
    // modal displacements matrix
    pub(crate) c: DMatrix<if64>,
    // static solution gain matrix
    pub(crate) g_ssol: Option<DMatrix<f64>>,
    // static gain mismatch compensation scheme
    pub(crate) static_gain_mismatch: Option<StaticGainCompensation>,
    // eigen frequencies
    pub(crate) w: Vec<f64>,
    // damping coefficient
    pub(crate) z: f64,
    // optical sensitivity matrix
    pub(crate) optical_senses: Option<DMatrix<f64>>,
}

/// FEM structural dynamic model builder
#[derive(Debug, Default)]
pub struct StructuralBuilder {
    // inputs: Vec<String>,
    // outputs: Vec<String>,
    // z: f64,
    pub(crate) built: Structural,
    pub(crate) min_eigen_frequency: Option<f64>,
    pub(crate) max_eigen_frequency: Option<f64>,
    pub(crate) file_name: String,
    // static_gain_mismatch: Option<StaticGainCompensation>,
}
impl StructuralBuilder {
    /// Sets the FEM modal damping coefficient
    pub fn damping(mut self, z: f64) -> Self {
        self.built.z = z;
        self
    }
    /// Truncates the eigen frequencies to and including `max_eigen_frequency`
    ///
    /// The number of modes is set accordingly
    pub fn max_eigen_frequency(mut self, max_eigen_frequency: Option<f64>) -> Self {
        self.max_eigen_frequency = max_eigen_frequency;
        self
    }
    /// Drops the eigen frequencies less than `min_eigen_frequency`
    ///
    /// The number of modes is set accordingly
    pub fn min_eigen_frequency(mut self, min_eigen_frequency: Option<f64>) -> Self {
        self.min_eigen_frequency = min_eigen_frequency;
        self
    }
    /// Sets the filename where [Structural] is seralize to
    pub fn filename<S: Into<String>>(mut self, file_name: S) -> Self {
        self.file_name = file_name.into();
        self
    }
    pub fn optical_sensitivities(mut self, mat: Option<DMatrix<f64>>) -> Self {
        self.built.optical_senses = mat;
        self
    }
    /* /// Enables the compensation of the static gain mismatch
    ///
    /// An optional delay `s``:w` may be added
    fn enable_static_gain_mismatch_compensation(mut self, maybe_delay: Option<f64>) -> Self {
        self.static_gain_mismatch = Some(Default::default());
        if let Some(value) = maybe_delay {
            self.static_gain_mismatch
                .as_mut()
                .and_then(|sgm| sgm.delay.replace(value));
        }
        self
    } */
    fn new(inputs: Vec<String>, outputs: Vec<String>) -> Self {
        let built = Structural {
            inputs,
            outputs,
            z: 2. / 100.,
            ..Default::default()
        };
        Self {
            built,
            file_name: "structural".into(),
            ..Default::default()
        }
    }
    /// Builds the [Structural] model
    pub fn build(self) -> Result<Structural> {
        // let repo = env::var("DATA_REPO").unwrap_or_else(|_| ".".to_string());
        // let path = Path::new(&repo).join(self.file_name).with_extension("bin");
        // if let Ok(file) = File::open(&path) {
        //     println!("loading structural from {:?}", path);
        //     let buffer = BufReader::new(file);
        //     let this: Structural =
        //         bincode::serde::decode_from_reader(buffer, bincode::config::standard())?;
        //     if !(this.inputs == self.inputs && this.outputs == self.outputs) {
        //         return Err(StructuralError::IOMismatch(
        //             path.to_str().unwrap().to_string(),
        //         ));
        //     }
        //     Ok(this)
        // } else {
        println!("building structural from FEM");
        let mut fem = FEM::from_env()?;
        println!("{fem}");

        fem.switch_inputs(Switch::Off, None)
            .switch_inputs_by_name(self.built.inputs.clone(), Switch::On)?
            .switch_outputs(Switch::Off, None)
            .switch_outputs_by_name(self.built.outputs.clone(), Switch::On)?;
        let b = DMatrix::<f64>::from_row_slice(fem.n_modes(), fem.n_inputs(), &fem.inputs2modes())
            .map(|x| Complex::new(x, 0f64));
        let c =
            DMatrix::<f64>::from_row_slice(fem.n_outputs(), fem.n_modes(), &fem.modes2outputs())
                .map(|x| Complex::new(x, 0f64));
        let g_ssol = fem.reduced_static_gain();
        let w = fem.eigen_frequencies_to_radians();

        // self.static_gain_mismatch.as_mut().map(|sgm| {
        //     let g_dsol = fem.static_gain();
        //     let delta_g = g_ssol.as_ref().expect("failed to get FEM static gain") - g_dsol;
        //     sgm.delta_gain = delta_g.map(|x| Complex::new(x, 0f64));
        // });

        let q = match (self.min_eigen_frequency, self.max_eigen_frequency) {
            (Some(min), Some(max)) => Some((
                fem.eigen_frequencies
                    .iter()
                    .copied()
                    .enumerate()
                    .find(|(_, f)| *f >= min)
                    .unwrap_or_default()
                    .0,
                fem.eigen_frequencies
                    .iter()
                    .copied()
                    .filter_map(|f| (f >= min && f <= max).then(|| f))
                    .enumerate()
                    .last()
                    .unwrap_or_default()
                    .0
                    + 1,
            )),
            (None, Some(max)) => Some((
                0,
                fem.eigen_frequencies
                    .iter()
                    .copied()
                    .filter_map(|f| (f <= max).then(|| f))
                    .enumerate()
                    .last()
                    .unwrap_or_default()
                    .0
                    + 1,
            )),
            (Some(min), None) => {
                let s = fem
                    .eigen_frequencies
                    .iter()
                    .copied()
                    .enumerate()
                    .find(|(_, f)| *f >= min)
                    .unwrap_or_default()
                    .0;
                Some((s, fem.eigen_frequencies.len() - s))
            }
            (None, None) => None,
        };

        Ok(if let Some((s, n)) = q {
            Structural {
                b: b.rows(s, n).into_owned(),
                c: c.columns(s, n).into_owned(),
                g_ssol,
                w: w[s..s + n].to_vec(),
                ..self.built
            }
        } else {
            Structural {
                b,
                c,
                g_ssol,
                w,
                ..self.built
            }
        })
        // let file = File::create(&path)?;
        // let mut buffer = BufWriter::new(file);
        // bincode::serde::encode_into_std_write(&this, &mut buffer, bincode::config::standard())?;
        // println!("structural save to {:?}", path);
        // Ok(this)
        // }
    }
}
impl Structural {
    /// Creates a [Structural] builder
    pub fn builder(inputs: Vec<String>, outputs: Vec<String>) -> StructuralBuilder {
        StructuralBuilder::new(inputs, outputs)
    }
    /// Returns a [view](https://docs.rs/nalgebra/latest/nalgebra/base/struct.Matrix.html#method.view) of the static gain
    pub fn static_gain(
        &self,
        ij: (usize, usize),
        nm: (usize, usize),
    ) -> Option<DMatrixView<'_, f64>> {
        self.g_ssol.as_ref().map(|g| g.view(ij, nm))
    }
    /// Returns the eigen frequencies in Hz
    pub fn eigen_frequencies_hz(&self) -> Vec<f64> {
        self.w
            .iter()
            .map(|x| *x * 0.5 * consts::FRAC_1_PI)
            .collect()
    }
}

impl Display for Structural {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "GMT structural dynamic model:")?;
        writeln!(f, " + inputs: {:?}", self.inputs)?;
        writeln!(f, " + outputs: {:?}", self.outputs)?;
        writeln!(
            f,
            " + eigen frequencies: ({:.3},{:.3})Hz",
            0.5 * self.w[0] * consts::FRAC_1_PI,
            0.5 * self.w.last().unwrap() * consts::FRAC_1_PI
        )?;
        writeln!(f, " + damping: {:}%", self.z * 1e2)?;
        writeln!(f, " + B matrix {:?}", self.b.shape())?;
        writeln!(f, " + C matrix {:?}", self.c.shape())?;
        if let Some(g) = self.g_ssol.as_ref() {
            writeln!(f, " + static gain matrix {:?}", g.shape())?;
        }
        Ok(())
    }
}

impl FrequencyResponse for Structural {
    type Output = DMatrix<Complex<f64>>;

    /// *Dynamics and Control of Structures, W.K. Gawronsky*, p.17-18, Eqs.(2.21)-(2.22)
    fn j_omega(&self, jw: if64) -> Self::Output {
        let zeros = DMatrix::<Complex<f64>>::zeros(self.c.nrows(), self.b.ncols());
        let fr = self
            .c
            .column_iter()
            .zip(self.b.row_iter())
            .zip(&self.w)
            .fold(zeros, |a, ((c, b), wi)| {
                let mut cb = c * b;
                let ode = wi * wi + jw * jw + 2f64 * self.z * wi * jw;
                cb /= ode;
                a + cb
            });
        let fr = match &self.static_gain_mismatch {
            Some(StaticGainCompensation {
                delay: None,
                delta_gain,
            }) => fr + delta_gain,
            Some(StaticGainCompensation {
                delay: Some(t_s),
                delta_gain,
            }) => fr + (delta_gain * (-jw * t_s).exp()),
            None => fr,
        };
        if let Some(mat) = self.optical_senses.as_ref() {
            mat.map(|x| Complex::from(x)) * fr
        } else {
            fr
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::frequency_response::Frequencies;

    use super::*;

    #[test]
    fn mount() {
        let structural = Structural::builder(
            vec!["OSS_ElDrive_Torque".to_string()],
            vec!["OSS_ElEncoder_Angle".to_string()],
        )
        .build()
        .unwrap();

        let tf = structural.frequency_response(1f64);
        println!("{}", tf[0]);
    }

    #[test]
    fn mount_el_tf() {
        let structural = Structural::builder(
            vec!["OSS_ElDrive_Torque".to_string()],
            vec!["OSS_ElEncoder_Angle".to_string()],
        )
        .build()
        .unwrap();

        let tf = structural.frequency_response(Frequencies::logspace(0.1, 100., 1000));
        println!("{}", tf[0]);

        // let mut file = File::create("mount_el_tf.pkl").unwrap();
        // serde_pickle::to_writer(&mut file, &(nu, tf), Default::default()).unwrap();
    }

    #[test]
    fn mount_el_tf_dc() {
        let structural = Structural::builder(
            vec!["OSS_ElDrive_Torque".to_string()],
            vec!["OSS_ElEncoder_Angle".to_string()],
        )
        // .enable_static_gain_mismatch_compensation(Some(1. / 8e3))
        .build()
        .unwrap();

        let tf = structural.frequency_response(Frequencies::logspace(0.1, 4e3, 1000));
        //println!("{:?}", nu);
        println!("{}", tf[0]);

        // let mut file = File::create("mount_el_tf_dc_full-sampling_delay.pkl").unwrap();
        // serde_pickle::to_writer(&mut file, &(nu, tf), Default::default()).unwrap();
    }

    #[test]
    fn mount_el_tf_linspace() {
        let structural = Structural::builder(
            vec!["OSS_ElDrive_Torque".to_string()],
            vec!["OSS_ElEncoder_Angle".to_string()],
        )
        .build()
        .unwrap();

        let tf = structural.frequency_response(Frequencies::LinSpace {
            lower: 1f64,
            upper: 10f64,
            n: 2,
        });
        // println!("{:?}", nu);
        println!("{}", tf[0]);

        // let sys = Sys::from((nu, tf));
        // dbg!(sys);
    }
}
