use clap::Parser;
use gmt_fem_frequency_response::{
    BuilderTrait, cli::Cli, data::TransferFunctionData, frequency_response::FrequencyResponse,
    structural::Structural,
};

fn main() -> anyhow::Result<()> {
    let args: Cli = Cli::parse();

    let model = Structural::builder(
        args.inputs.iter().map(|io| io.name()).collect(),
        args.outputs.iter().map(|io| io.name()).collect(),
    )
    .damping(args.structural_damping)
    .build()?;

    let nu = args.frequencies.clone();
    let frequency_response = model.frequency_response(nu);

    TransferFunctionData::from(&args)
        .add_response(frequency_response)
        .dump(args.filename)?;

    Ok(())
}
