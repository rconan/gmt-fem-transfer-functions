use std::time::Instant;

use clap::Parser;
use gmt_fem_frequency_response::{
    Cli,
    data::{Extremum, TransferFunctionData},
    structural::{Structural, StructuralFrequencyResponse},
};

fn main() -> anyhow::Result<()> {
    let args: Cli = Cli::parse();

    let model = Structural::try_from(&args)?;
    println!("{model}");

    let now = Instant::now();
    let frequency_response = if args.svd {
        println!("computing frequency response SVD");
        model.frequency_response_svd(&args.frequencies)
    } else {
        model.frequency_response(&args.frequencies)
    };
    println!(
        "frequency response computed in {:.3}s",
        now.elapsed().as_secs_f64()
    );
    println!("{frequency_response}");

    let mut ex = frequency_response.extrema(None);
    ex.sort_by(|Extremum { y: a, .. }, Extremum { y: b, .. }| b.partial_cmp(a).unwrap());
    println!("Frequency response extrema:");
    ex.iter()
        .take(5)
        .enumerate()
        .for_each(|(i, Extremum { x, y })| println!(" {:2}: {:8.2} {:.3e}", i + 1, x, y));

    TransferFunctionData::from(&args)
        .add_structural(&model)
        .add_response(frequency_response)
        .dump(args.filename)?;

    Ok(())
}
