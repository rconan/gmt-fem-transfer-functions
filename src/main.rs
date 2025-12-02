use std::time::Instant;

use clap::Parser;
use gmt_fem_frequency_response::{
    Cli, data::TransferFunctionData, frequency_response::FrequencyResponse, structural::Structural,
};

fn main() -> anyhow::Result<()> {
    let args: Cli = Cli::parse();

    let model = Structural::try_from(&args)?;
    println!("{model}");

    let now = Instant::now();
    let frequency_response = model.frequency_response(&args.frequencies);
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
