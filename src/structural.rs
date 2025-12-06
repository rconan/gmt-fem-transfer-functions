//! FEM structural dynamic model

use std::{f64::consts, fmt::Display};

#[cfg(feature = "faer")]
use faer::{Mat, MatRef};
use gmt_dos_clients_fem::{Model, Switch};
use gmt_fem::FEM;
#[cfg(feature = "nalgebra")]
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
    #[cfg(feature = "nalgebra")]
    pub(crate) delta_gain: DMatrix<if64>,
    #[cfg(feature = "faer")]
    pub(crate) delta_gain: Mat<if64>,
}
impl Default for StaticGainCompensation {
    fn default() -> Self {
        Self {
            delay: Default::default(),
            #[cfg(feature = "nalgebra")]
            delta_gain: DMatrix::<if64>::zeros(1, 1),
            #[cfg(feature = "faer")]
            delta_gain: Mat::<if64>::zeros(1, 1),
        }
    }
}

/// FEM structural dynamic model
#[derive(Debug, Deserialize, Serialize)]
pub struct Structural {
    // inputs labels
    pub(crate) inputs: Vec<String>,
    // outputs labels
    pub(crate) outputs: Vec<String>,
    // modal forces matrix
    #[cfg(feature = "nalgebra")]
    pub(crate) b: DMatrix<f64>,
    #[cfg(feature = "faer")]
    pub(crate) b: Mat<if64>,
    // modal displacements matrix
    #[cfg(feature = "nalgebra")]
    pub(crate) c: DMatrix<f64>,
    #[cfg(feature = "faer")]
    pub(crate) c: Mat<if64>,
    // static solution gain matrix
    #[cfg(feature = "nalgebra")]
    pub(crate) g_ssol: Option<DMatrix<f64>>,
    #[cfg(feature = "faer")]
    pub(crate) g_ssol: Option<Mat<f64>>,
    // static gain mismatch compensation scheme
    pub(crate) static_gain_mismatch: Option<StaticGainCompensation>,
    // eigen frequencies
    pub(crate) w: Vec<f64>,
    // damping coefficient
    pub(crate) z: f64,
    // optical sensitivity matrix
    #[cfg(feature = "nalgebra")]
    pub(crate) optical_senses: Option<DMatrix<f64>>,
    #[cfg(feature = "faer")]
    pub(crate) optical_senses: Option<Mat<Complex<f64>>>,
}

#[cfg(feature = "nalgebra")]
impl Default for Structural {
    fn default() -> Self {
        Self {
            inputs: Default::default(),
            outputs: Default::default(),
            b: Default::default(),
            c: Default::default(),
            g_ssol: Default::default(),
            static_gain_mismatch: Default::default(),
            w: Default::default(),
            z: Default::default(),
            optical_senses: Default::default(),
        }
    }
}
#[cfg(feature = "faer")]
impl Default for Structural {
    fn default() -> Self {
        Self {
            inputs: Default::default(),
            outputs: Default::default(),
            b: Mat::new(),
            c: Mat::new(),
            g_ssol: Default::default(),
            static_gain_mismatch: Default::default(),
            w: Default::default(),
            z: Default::default(),
            optical_senses: Default::default(),
        }
    }
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
    /// Sets the optical sensitivity matrix
    #[cfg(feature = "nalgebra")]
    pub fn optical_sensitivities(mut self, mat: Option<DMatrix<f64>>) -> Self {
        self.built.optical_senses = mat;
        self
    }
    #[cfg(feature = "faer")]
    pub fn optical_sensitivities(mut self, mat: Option<Mat<Complex<f64>>>) -> Self {
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
        #[cfg(feature = "nalgebra")]
        let b = DMatrix::<f64>::from_row_slice(fem.n_modes(), fem.n_inputs(), &fem.inputs2modes());
        #[cfg(feature = "faer")]
        let b = MatRef::<Complex<f64>>::from_row_major_slice(
            &fem.inputs2modes()
                .into_iter()
                .map(|x| Complex::new(x, 0f64))
                .collect::<Vec<_>>(),
            fem.n_modes(),
            fem.n_inputs(),
        )
        .to_owned();
        #[cfg(feature = "nalgebra")]
        let c =
            DMatrix::<f64>::from_row_slice(fem.n_outputs(), fem.n_modes(), &fem.modes2outputs());
        #[cfg(feature = "faer")]
        let c = MatRef::<Complex<f64>>::from_row_major_slice(
            &fem.modes2outputs()
                .into_iter()
                .map(|x| Complex::new(x, 0f64))
                .collect::<Vec<_>>(),
            fem.n_outputs(),
            fem.n_modes(),
        )
        .to_owned();
        #[cfg(feature = "nalgebra")]
        let g_ssol = fem.reduced_static_gain();
        #[cfg(feature = "faer")]
        let g_ssol = None;
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
                #[cfg(feature = "nalgebra")]
                b: b.rows(s, n).into_owned(),
                #[cfg(feature = "faer")]
                b: b.subrows(s, n).to_owned(),
                #[cfg(feature = "nalgebra")]
                c: c.columns(s, n).into_owned(),
                #[cfg(feature = "faer")]
                c: c.subcols(s, n).to_owned(),
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
    #[cfg(feature = "nalgebra")]
    pub fn static_gain(
        &self,
        ij: (usize, usize),
        nm: (usize, usize),
    ) -> Option<DMatrixView<'_, f64>> {
        self.g_ssol.as_ref().map(|g| g.view(ij, nm))
    }
    #[cfg(feature = "faer")]
    pub fn static_gain(&self, _ij: (usize, usize), _nm: (usize, usize)) -> Option<MatRef<'_, f64>> {
        None
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

#[cfg(feature = "nalgebra")]
impl FrequencyResponse for Structural {
    type Output = DMatrix<Complex<f64>>;

    /// *Dynamics and Control of Structures, W.K. Gawronsky*, p.17-18, Eqs.(2.21)-(2.22)
    fn j_omega(&self, jw: if64) -> Self::Output {
        let zeros = DMatrix::<Complex<f64>>::zeros(self.c.nrows(), self.b.ncols());
        let mut cb = DMatrix::<f64>::zeros(self.c.nrows(), self.b.ncols());
        let mut ccb = DMatrix::<if64>::zeros(self.c.nrows(), self.b.ncols());
        let fr = self
            .c
            .column_iter()
            .zip(self.b.row_iter())
            .zip(&self.w)
            .fold(zeros, |a, ((c, b), wi)| {
                let ode = 1f64 / (wi * wi + jw * jw + 2f64 * self.z * wi * jw);
                // let now = std::time::Instant::now();
                // let cb = (c * b);
                c.mul_to(&b, &mut cb);
                // cb_rt += now.elapsed().as_micros();
                // cb /= ode;
                ccb.zip_apply(&cb, |l, r| *l = Complex::from(r) * ode);
                a + &ccb //.map(|x| Complex::from(x) * ode)
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
#[cfg(feature = "faer")]
impl FrequencyResponse for Structural {
    type Output = Mat<Complex<f64>>;

    /// *Dynamics and Control of Structures, W.K. Gawronsky*, p.17-18, Eqs.(2.21)-(2.22)
    fn j_omega(&self, jw: if64) -> Self::Output {
        use faer::{Accum, diag::DiagRef, get_global_parallelism, linalg::matmul::matmul};
        let mut fr = Mat::<Complex<f64>>::zeros(self.c.nrows(), self.b.ncols());
        let rode: Vec<_> = self
            .w
            .iter()
            .map(|wi| wi * wi + jw * jw + 2f64 * self.z * wi * jw)
            .map(|ode| 1f64 / ode)
            .collect();
        let d = DiagRef::from_slice(&rode);
        matmul(
            &mut fr,
            Accum::Replace,
            &self.c,
            d * &self.b,
            1f64.into(),
            get_global_parallelism(),
        );
        // let fr = match &self.static_gain_mismatch {
        //     Some(StaticGainCompensation {
        //         delay: None,
        //         delta_gain,
        //     }) => fr + delta_gain,
        //     Some(StaticGainCompensation {
        //         delay: Some(t_s),
        //         delta_gain,
        //     }) => fr + (delta_gain * (-jw * t_s).exp()),
        //     None => fr,
        // };
        if let Some(mat) = self.optical_senses.as_ref() {
            mat * fr
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
