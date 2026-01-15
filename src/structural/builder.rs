#[cfg(feature = "faer")]
use faer::{Mat, MatRef};
#[cfg(feature = "nalgebra")]
use nalgebra::DMatrix;
use gmt_dos_clients_fem::{Model, Switch};
use gmt_fem::FEM;

#[cfg(feature = "faer")]
use crate::if64;

pub use super::Structural;

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
    pub fn optical_sensitivities(mut self, mat: Option<Mat<if64>>) -> Self {
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
    pub fn new(inputs: Vec<String>, outputs: Vec<String>) -> Self {
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
    pub fn build(self) -> super::Result<Structural> {
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
        let b = DMatrix::<f64>::from_row_slice(
            fem.n_modes(),
            fem.n_inputs(),
            &fem.named_inputs_to_modes(&self.built.inputs)?.unwrap(),
        );
        #[cfg(feature = "faer")]
        let b = MatRef::<if64>::from_row_major_slice(
            &fem.named_inputs_to_modes(&self.built.inputs)?
                .unwrap()
                .into_iter()
                .map(|x| if64::new(x, 0f64))
                .collect::<Vec<_>>(),
            fem.n_modes(),
            fem.n_inputs(),
        )
        .to_owned();
        #[cfg(feature = "nalgebra")]
        let c = DMatrix::<f64>::from_row_slice(
            fem.n_outputs(),
            fem.n_modes(),
            &fem.modes_to_named_outputs(&self.built.outputs)?.unwrap(),
        );
        #[cfg(feature = "faer")]
        let c = MatRef::<if64>::from_row_major_slice(
            &fem.modes_to_named_outputs(&self.built.outputs)?
                .unwrap()
                .into_iter()
                .map(|x| if64::new(x, 0f64))
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
        //     sgm.delta_gain = delta_g.map(|x| if64::new(x, 0f64));
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
