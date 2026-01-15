use std::time::Instant;

use clap::Parser;
use gmt_fem_frequency_response::{
    Cli,
    data::TransferFunctionData,
    if64,
    structural::{Structural, StructuralFrequencyResponse},
};
#[cfg(feature = "faer")]
use gmt_fem_frequency_response::{data::Extremum, structural::StructuralFrequencyResponseSvd};

fn main() -> anyhow::Result<()> {
    #[allow(unused_mut)]
    let mut args: Cli = Cli::parse();

    let model = Structural::try_from(&args)?;
    println!("{model}");

    let now = Instant::now();
    if args.svd {
        #[cfg(feature = "nalgebra")]
        {
            return Err(anyhow::anyhow!(
                "SVD requires the faer feature, it is not implemented for the nalgebra feature"
                    .to_string(),
            ));
        }
        #[cfg(feature = "faer")]
        {
            let (u, v) = match args
                .uv
                .take()
                .map(|svd| svd.to_lowercase())
                .as_ref()
                .map(|x| x.as_str())
            {
                Some("uv") => (true, true),
                Some("u") => (true, false),
                Some("v") => (false, true),
                Some(other) => panic!(r#"found svd argument {other}, expected "u", "v" or "uv"#),
                None => (false, false),
            };
            println!("computing frequency response SVD");
            let frequency_response = model.frequency_response_svd(&args.frequencies, u, v);
            println!(
                "frequency response computed in {:.3}s",
                now.elapsed().as_secs_f64()
            );
            println!("{frequency_response}");

            let mut ex = frequency_response.extrema(None);
            ex.sort_by(|Extremum { y: a, .. }, Extremum { y: b, .. }| b.partial_cmp(a).unwrap());
            println!("Sorted frequency response extrema:");
            ex.iter()
                .take(10)
                .for_each(|Extremum { i, x, y }| println!(" {:5}: {:8.2} {:.3e}", i, x, y));

            TransferFunctionData::<f64>::from(&args)
                .add_structural(&model, args.b, args.c)?
                .add_response(frequency_response)
                .dump(args.filename)?;
        }
    } else {
        let frequency_response = model.frequency_response(&args.frequencies);
        println!(
            "frequency response computed in {:.3}s",
            now.elapsed().as_secs_f64()
        );
        println!("{frequency_response}");

        TransferFunctionData::<if64>::from(&args)
            .add_structural(&model, args.b, args.c)?
            .add_response(frequency_response)
            .dump(args.filename)?;
    };

    Ok(())
}
