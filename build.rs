use std::{env, fs, path::Path};

fn main() {
    let _ = gmt_fem_code_builder::rustc_config(env!("CARGO_PKG_NAME"), None);
    match gmt_fem_code_builder::io_names(env!("CARGO_PKG_NAME")) {
        Ok((inputs, outputs)) => {
            let inputs_variants: Vec<_> = inputs.iter().map(|io| io.variant()).collect();
            let inputs_variant_to_name: Vec<_> = inputs
                .iter()
                .map(|io| format!(r#"Self::{} => "{}""#, io.variant(), io.name))
                .collect();
            let outputs_variants: Vec<_> = outputs.iter().map(|io| io.variant()).collect();
            let outputs_variant_to_name: Vec<_> = outputs
                .iter()
                .map(|io| format!(r#"Self::{} => "{}""#, io.variant(), io.name))
                .collect();
            let out_dir = env::var("OUT_DIR").unwrap();
            let path = Path::new(&out_dir).join("fem_io.rs");
            match fs::write(
                &path,
                format!(
                    r#"
#[doc = "FEM inputs"]
#[derive(Clone, ValueEnum)]
pub enum Inputs {{
    {}
}}
impl Inputs {{
    pub fn name(&self) -> String {{
        match self {{
            {}
        }}.to_string()
    }}
}}
#[doc = "FEM outputs"]
#[derive(Clone, ValueEnum)]
pub enum Outputs {{
    {},
    TipTilt,
    SegmentTipTilt,
    SegmentPiston,
}}
impl Outputs {{
    pub fn name(&self) -> String {{
        match self {{
            {},
            Self::TipTilt => "tip-tilt",
            Self::SegmentTipTilt => "segment_tip-tilt",
            Self::SegmentPiston => "segment_piston",
        }}.to_string()
    }}
}}
                "#,
                    inputs_variants.join(",\n    "),
                    inputs_variant_to_name.join(",\n            "),
                    outputs_variants.join(",\n    "),
                    outputs_variant_to_name.join(",\n            "),
                ),
            ) {
                Ok(_) => {
                    println!("cargo::warning=FEM IO written to {:?}", &path);
                }
                Err(e) => println!("cargo::error={e}"),
            };
        }
        Err(e) => println!("cargo::error={e}"),
    };
}
