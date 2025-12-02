use std::time::Instant;

use clap::Parser;
use gmt_fem_frequency_response::{
    Cli, data::TransferFunctionData, frequency_response::FrequencyResponse, structural::Structural,
};

fn main() -> anyhow::Result<()> {
    let args: Cli = Cli::parse();

    let model = Structural::builder(args.fem_inputs(), args.fem_outputs())
        .damping(args.structural_damping)
        .min_eigen_frequency(args.eigen_frequency_min)
        .max_eigen_frequency(args.eigen_frequency_max)
        .build()?;
    println!("{model}");

    let nu = args.frequencies.clone();
    let now = Instant::now();
    let frequency_response = model.frequency_response(nu);
    println!(
        "frequency response computed in {:.3}s",
        now.elapsed().as_secs_f64()
    );
    println!("{frequency_response}");

    TransferFunctionData::from(&args)
        .add_response(frequency_response)
        .dump(args.filename)?;

    Ok(())
}
