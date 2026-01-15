//! FEM structural dynamic model
use std::{f64::consts, fmt::Display};

mod builder;
pub use builder::StructuralBuilder;

#[cfg(feature = "faer")]
use faer::{Mat, MatRef};
use gmt_dos_clients_fem::{Model, Switch};
use gmt_fem::FEM;
#[cfg(feature = "nalgebra")]
use nalgebra::{DMatrix, DMatrixView};
use serde::{Deserialize, Serialize};

#[cfg(feature = "faer")]
use crate::frequency_response::{FrequencyResponseSvd, JOmegaSvd};
use crate::{
    frequency_response::{Frequencies, FrequencyResponse, JOmega},
    if64,
};

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
    pub(crate) optical_senses: Option<Mat<if64>>,
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
    /// Return the model natural frequencies
    pub fn natural_frequencies(
        &self,
        min: Option<f64>,
        max: Option<f64>,
        index: Option<Vec<usize>>,
    ) -> Frequencies {
        Frequencies::Set {
            values: {
                if let Some(index) = index {
                    index
                        .into_iter()
                        .map(|i| self.w[i])
                        .map(|x| 0.5 * x * std::f64::consts::FRAC_1_PI)
                        .collect()
                } else {
                    self.w
                        .iter()
                        .copied()
                        .map(|x| 0.5 * x * std::f64::consts::FRAC_1_PI)
                        .filter_map(|x| {
                            if let Some(min) = min
                                && x < min
                            {
                                None
                            } else {
                                Some(x)
                            }
                        })
                        .filter_map(|x| {
                            if let Some(max) = max
                                && x > max
                            {
                                None
                            } else {
                                Some(x)
                            }
                        })
                        .collect()
                }
            },
        }
    }
    /// Returns the inputs nodes \[x,y,z\]
    pub fn inputs_nodes(&self) -> Result<Vec<f64>> {
        let mut fem = FEM::from_env()?;
        fem.switch_inputs(Switch::Off, None)
            .switch_inputs(Switch::Off, None)
            .switch_inputs_by_name(self.inputs.clone(), Switch::On)?;
        Ok(self
            .inputs
            .iter()
            .flat_map(|input| {
                let get_out =
                    Box::<dyn gmt_dos_clients_fem::fem_io::GetIn>::try_from(input.clone()).unwrap();
                let idx = get_out.position(&fem.inputs).unwrap();
                fem.inputs[idx]
                    .as_ref()
                    .map(|i| i.get_by(|i| i.properties.location.clone()))
                    .unwrap()
            })
            .flatten()
            .collect())
    }
    /// Returns the outputs nodes \[x,y,z\]
    pub fn outputs_nodes(&self) -> Result<Vec<f64>> {
        let mut fem = FEM::from_env()?;
        fem.switch_inputs(Switch::Off, None)
            .switch_outputs(Switch::Off, None)
            .switch_outputs_by_name(self.outputs.clone(), Switch::On)?;
        Ok(self
            .outputs
            .iter()
            .flat_map(|output| {
                let get_out =
                    Box::<dyn gmt_dos_clients_fem::fem_io::GetOut>::try_from(output.clone())
                        .unwrap();
                let idx = get_out.position(&fem.outputs).unwrap();
                fem.outputs[idx]
                    .as_ref()
                    .map(|i| i.get_by(|i| i.properties.location.clone()))
                    .unwrap()
            })
            .flatten()
            .collect())
    }
    /// Returns the force to mode matrix
    #[cfg(feature = "faer")]
    pub fn force_to_mode(&self) -> Mat<f64> {
        let mut iter = self.b.col_iter().flat_map(|c| c.iter().map(|x| x.re));
        Mat::<f64>::from_fn(self.b.nrows(), self.b.ncols(), |_, _| iter.next().unwrap())
    }
    #[cfg(feature = "nalgebra")]
    pub fn force_to_mode(&self) -> DMatrix<f64> {
        self.b.clone()
    }
    /// Returns the mode to displacement matrix
    #[cfg(feature = "faer")]
    pub fn mode_to_displacement(&self) -> Mat<f64> {
        let mut iter = self.c.col_iter().flat_map(|c| c.iter().map(|x| x.re));
        Mat::<f64>::from_fn(self.c.nrows(), self.c.ncols(), |_, _| iter.next().unwrap())
    }
    #[cfg(feature = "nalgebra")]
    pub fn mode_to_displacement(&self) -> DMatrix<f64> {
        self.c.clone()
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
impl JOmega for Structural {
    type Output = DMatrix<if64>;

    /// *Dynamics and Control of Structures, W.K. Gawronsky*, p.17-18, Eqs.(2.21)-(2.22)
    fn j_omega(&self, jw: if64) -> Self::Output {
        let zeros = DMatrix::<if64>::zeros(self.c.nrows(), self.b.ncols());
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
                ccb.zip_apply(&cb, |l, r| *l = if64::from(r) * ode);
                a + &ccb //.map(|x| if64::from(x) * ode)
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
            mat.map(|x| if64::from(x)) * fr
        } else {
            fr
        }
    }
}
#[cfg(feature = "faer")]
impl JOmega for Structural {
    type Output = Mat<if64>;

    /// *Dynamics and Control of Structures, W.K. Gawronsky*, p.17-18, Eqs.(2.21)-(2.22)
    fn j_omega(&self, jw: if64) -> Self::Output {
        use faer::{Accum, diag::DiagRef, get_global_parallelism, linalg::matmul::matmul};
        let mut fr = Mat::<if64>::zeros(self.c.nrows(), self.b.ncols());
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
#[cfg(feature = "faer")]
impl JOmegaSvd for Structural {
    type Svd = Mat<f64>;
    fn j_omega_svd(
        &self,
        jw: if64,
        u: bool,
        v: bool,
    ) -> (Self::Svd, Option<Self::Svd>, Option<Self::Svd>) {
        use crate::data::Cartesian2Polar;

        let mat = self.j_omega(jw).magnitude();
        let svd = mat.svd().unwrap();
        let s = svd.S().column_vector().as_mat().to_owned();
        match (u, v) {
            (true, true) => (s, Some(svd.U().to_owned()), Some(svd.V().to_owned())),
            (true, false) => (s, Some(svd.U().to_owned()), None),
            (false, true) => (s, None, Some(svd.V().to_owned())),
            (false, false) => (s, None, None),
        }
    }
}

pub trait StructuralFrequencyResponse: FrequencyResponse {
    fn frequency_response<T: Into<crate::frequency_response::Frequencies>>(
        &self,
        nu: T,
    ) -> crate::data::FrequencyResponseVec<Self::Output>
    where
        <Self as JOmega>::Output: crate::data::Cartesian2Polar + Send,
        <<Self as JOmega>::Output as crate::data::Cartesian2Polar>::Output: Send,
        Self: Sync;
}
#[cfg(feature = "faer")]
pub trait StructuralFrequencyResponseSvd: FrequencyResponseSvd {
    fn frequency_response_svd<T: Into<crate::frequency_response::Frequencies>>(
        &self,
        nu: T,
        u: bool,
        v: bool,
    ) -> crate::data::FrequencyResponseVec<Self::Svd>
    where
        <Self as JOmegaSvd>::Svd: crate::data::Cartesian2Polar + Send,
        <<Self as JOmegaSvd>::Svd as crate::data::Cartesian2Polar>::Output: Send,
        Self: Sync;
}

impl StructuralFrequencyResponse for Structural {
    fn frequency_response<T: Into<crate::frequency_response::Frequencies>>(
        &self,
        nu: T,
    ) -> crate::data::FrequencyResponseVec<Self::Output>
    where
        <Self as JOmega>::Output: crate::data::Cartesian2Polar + Send,
        <<Self as JOmega>::Output as crate::data::Cartesian2Polar>::Output: Send,
        Self: Sync,
    {
        let frequencies: Frequencies = nu.into();
        let frequencies = if let Frequencies::Structural { min, max, index } = frequencies {
            self.natural_frequencies(min, max, index)
        } else {
            frequencies
        };
        <Self as FrequencyResponse>::frequency_response(&self, frequencies)
    }
}
#[cfg(feature = "faer")]
impl StructuralFrequencyResponseSvd for Structural {
    fn frequency_response_svd<T: Into<crate::frequency_response::Frequencies>>(
        &self,
        nu: T,
        u: bool,
        v: bool,
    ) -> crate::data::FrequencyResponseVec<Self::Svd>
    where
        <Self as JOmegaSvd>::Svd: crate::data::Cartesian2Polar + Send,
        <<Self as JOmegaSvd>::Svd as crate::data::Cartesian2Polar>::Output: Send,
        Self: Sync,
    {
        let frequencies: Frequencies = nu.into();
        let frequencies = if let Frequencies::Structural { min, max, index } = frequencies {
            self.natural_frequencies(min, max, index)
        } else {
            frequencies
        };
        <Self as FrequencyResponseSvd>::frequency_response_svd(&self, frequencies, u, v)
    }
}

#[cfg(test)]
mod tests {
    use crate::frequency_response::Frequencies;

    use super::{Structural, StructuralFrequencyResponse};

    #[test]
    fn mount() {
        let structural = Structural::builder(
            vec!["OSS_ElDrive_Torque".to_string()],
            vec!["OSS_ElEncoder_Angle".to_string()],
        )
        .build()
        .unwrap();

        let tf = structural.frequency_response(1f64);
        println!("{:?}", tf[0]);
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
        println!("{:?}", tf[0]);

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
        println!("{:?}", tf[0]);

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
        println!("{:?}", tf[0]);

        // let sys = Sys::from((nu, tf));
        // dbg!(sys);
    }
}
