use clap::Parser;
use gmt_fem_transfer_functions::{
    BuilderTrait, Inputs, Outputs,
    data::TransferFunctionData,
    frequency_response::{Frequencies, FrequencyResponse},
    structural::Structural,
};

#[derive(Parser)]
pub struct Cli {
    /// FEM inputs
    #[arg(short, long)]
    inputs: Vec<Inputs>,
    /// FEM outputs
    #[arg(short, long)]
    outputs: Vec<Outputs>,
    /// FEM modal damping coeffcient
    #[arg(short = 'z', long, default_value_t = 0.02f64)]
    structural_damping: f64,
    /// Frequencies [Hz]
    #[command(subcommand)]
    frequencies: Frequencies,
}

fn main() -> anyhow::Result<()> {
    let args: Cli = Cli::parse();

    let model = Structural::builder(
        args.inputs.iter().map(|io| io.name()).collect(),
        args.outputs.iter().map(|io| io.name()).collect(),
    )
    .damping(args.structural_damping)
    .build()?;

    // let nu = Frequencies::logspace(1f64, 250f64, 1000);
    let nu = Frequencies::from(1f64);
    let frequency_response = model.frequency_response(nu);

    let data = TransferFunctionData { frequency_response };
    data.dump("test.pkl")?;

    Ok(())
}
