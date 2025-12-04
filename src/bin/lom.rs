use gmt_fem_frequency_response::cli::Lom;
use gmt_lom::{Loader, LoaderTrait, OpticalSensitivities, OpticalSensitivity};
use std::{fs::File, io::Write};

fn main() -> anyhow::Result<()> {
    let sensitivities = Loader::<OpticalSensitivities>::default().load()?;

    let mut lom = Lom::default();

    vec![
        &sensitivities[OpticalSensitivity::TipTilt(vec![])],
        &sensitivities[OpticalSensitivity::SegmentTipTilt(vec![])],
        &sensitivities[OpticalSensitivity::SegmentPiston(vec![])],
    ]
    .iter()
    .for_each(|sens| {
        match sens {
            OpticalSensitivity::TipTilt(items) => {
                lom.tip_tilt = items.to_vec();
            }
            OpticalSensitivity::SegmentTipTilt(items) => {
                lom.segment_tip_tilt = items.to_vec();
            }
            OpticalSensitivity::SegmentPiston(items) => {
                lom.segment_piston = items.to_vec();
            }
            _ => {}
        };
    });

    let encoded = bincode::serde::encode_to_vec(&lom, bincode::config::standard())?;
    let compressed = lz4_flex::compress_prepend_size(&encoded);
    let mut file = File::create("src/lom.lz4").unwrap();
    file.write_all(&compressed).unwrap();

    Ok(())
}
