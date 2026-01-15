use crate::frequency_response::JOmega;
use crate::if64;

use std::f64::consts::PI;

const DPI: f64 = 2f64 * PI;

/// First order low-pass
///
/// *GMT-DOC-XXXX: ASM segment modal tranfer function*, Eq.(1)
#[derive(Debug)]
pub struct FirstOrderLowPass {
    corner_frequency_hz: f64,
}
impl FirstOrderLowPass {
    pub fn new() -> Self {
        Self {
            corner_frequency_hz: 4e3,
        }
    }
}
impl JOmega for FirstOrderLowPass {
    type Output = if64;
    fn j_omega(&self, jw: if64) -> Self::Output {
        jw / (1f64 + jw / (DPI * self.corner_frequency_hz))
    }
}

/// 4th-order bessel filter
///
/// *GMT-DOC-XXXX: ASM segment modal tranfer function*, Eq.(2)
#[derive(Debug)]
pub struct BesselFilter {
    w_bf: f64,
    beta: [f64; 5],
}
impl BesselFilter {
    pub fn new() -> Self {
        Self {
            w_bf: DPI * 2.2e3,
            beta: [1f64, 3.20108587, 4.39155033, 3.12393994, 1f64],
        }
    }
}
impl JOmega for BesselFilter {
    type Output = if64;
    fn j_omega(&self, jw: if64) -> Self::Output {
        let num = self.beta[0] * self.w_bf.powi(4);
        let denom = self
            .beta
            .iter()
            .enumerate()
            .fold(if64::new(0f64, 0f64), |a, (i, b)| {
                a + b * self.w_bf.powi(4 - i as i32) * jw.powi(i as i32)
            });
        num / denom
    }
}

/// Proportional-integral compensator
///
/// *GMT-DOC-XXXX: ASM segment modal tranfer function*, Eq.(3)
#[derive(Debug)]
pub struct PICompensator {
    kp: f64,
    ki: f64,
}
impl PICompensator {
    pub fn new() -> Self {
        Self { kp: 7e4, ki: 5e5 }
    }
}
impl JOmega for PICompensator {
    type Output = if64;
    fn j_omega(&self, jw: if64) -> Self::Output {
        self.kp + self.ki / jw
    }
}
#[cfg(test)]
mod tests {
    // use std::fs::File;

    use crate::frequency_response::{Frequencies, FrequencyResponse};

    use super::*;

    #[test]
    fn folp_tf() {
        let folp = FirstOrderLowPass::new();

        let tf = folp.frequency_response(Frequencies::logspace(1., 8e3, 1000));

        // let mut file = File::create("folp_tf.pkl").unwrap();
        // serde_pickle::to_writer(&mut file, &(nu, tf), Default::default()).unwrap();
    }

    #[test]
    fn bessel_tf() {
        let bessel = BesselFilter::new();

        let tf = bessel.frequency_response(Frequencies::logspace(1., 8e3, 1000));

        // let mut file = File::create("bessel_tf.pkl").unwrap();
        // serde_pickle::to_writer(&mut file, &(nu, tf), Default::default()).unwrap();
    }

    #[test]
    fn pic_tf() {
        let pic = PICompensator::new();

        let tf = pic.frequency_response(Frequencies::logspace(1., 8e3, 1000));

        // let mut file = File::create("pic_tf.pkl").unwrap();
        // serde_pickle::to_writer(&mut file, &(nu, tf), Default::default()).unwrap();
    }
}
