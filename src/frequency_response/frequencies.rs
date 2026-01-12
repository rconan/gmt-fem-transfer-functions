/// Frequency sampling options
///
/// The frequencies units is Hz
#[derive(Debug, Clone, clap::Subcommand)]
#[command(
    subcommand_help_heading = "Transfer functions sampling frequencies [Hz]",
    subcommand_value_name = "SAMPLING FREQUENCIES"
)]
pub enum Frequencies {
    /// a single frequency
    Single { value: f64 },
    /// logarithmic (log base 10) sampling of the interval `[lower,upper]` with `n` samples
    LogSpace {
        #[arg(short, long)]
        lower: f64,
        #[arg(short, long)]
        upper: f64,
        #[arg(short)]
        n: usize,
    },
    /// regular sampling of the interval `[lower,upper]` with `n` samples
    LinSpace {
        #[arg(short, long)]
        lower: f64,
        #[arg(short, long)]
        upper: f64,
        #[arg(short)]
        n: usize,
    },
    /// a given set of frequencies
    Set {
        #[arg(short, long)]
        values: Vec<f64>,
    },
    /// structural model natural frequencies
    Structural {
        /// minimum natural frequency [Hz]
        #[arg(long)]
        min: Option<f64>,
        /// maximum natural frequency [Hz]
        #[arg(long)]
        max: Option<f64>,
        /// natural frequency index
        #[arg(short, long)]
        index: Option<Vec<usize>>,
    },
}
impl From<f64> for Frequencies {
    fn from(value: f64) -> Self {
        Frequencies::Single { value }
    }
}
impl From<Vec<f64>> for Frequencies {
    fn from(values: Vec<f64>) -> Self {
        Frequencies::Set { values }
    }
}
impl From<&Vec<f64>> for Frequencies {
    fn from(values: &Vec<f64>) -> Self {
        Frequencies::Set {
            values: values.clone(),
        }
    }
}
impl From<&Self> for Frequencies {
    fn from(value: &Self) -> Self {
        value.clone()
    }
}
impl Frequencies {
    pub fn logspace(lower: f64, upper: f64, n: usize) -> Self {
        Self::LogSpace { lower, upper, n }
    }
    pub fn linspace(lower: f64, upper: f64, n: usize) -> Self {
        Self::LinSpace { lower, upper, n }
    }
}
