//! Frequency response data products

#[cfg(feature = "faer")]
use faer::Mat;
#[cfg(feature = "nalgebra")]
use nalgebra::{ComplexField, DMatrix};
use serde::Serialize;
use std::f64;

use crate::if64;

mod frequency;
mod transfer_function;
pub use frequency::{Extremum, FrequencyResponseData, FrequencyResponseVec};
pub use transfer_function::{ModalMatrix, TransferFunctionData, TransferFunctionDataError};

/// Matrix and scale size interface
pub trait Dims {
    type D: std::fmt::Debug + Serialize;
    fn size(&self) -> Self::D;
}

#[cfg(feature = "faer")]
impl Dims for Mat<f64> {
    type D = (usize, usize);

    fn size(&self) -> Self::D {
        self.shape()
    }
}
#[cfg(feature = "nalgebra")]
impl Dims for DMatrix<f64> {
    type D = (usize, usize);

    fn size(&self) -> Self::D {
        self.shape()
    }
}

impl Dims for f64 {
    type D = usize;

    fn size(&self) -> Self::D {
        1
    }
}

/// Cartesian to polar transformation interface
pub trait Cartesian2Polar {
    type Output: Dims + std::fmt::Debug + Serialize;
    fn magnitude(&self) -> Self::Output;
    fn phase(&self) -> Option<Self::Output>;
}

#[cfg(feature = "faer")]
impl Cartesian2Polar for Mat<if64> {
    type Output = Mat<f64>;

    fn magnitude(&self) -> Self::Output {
        let mut col_wise_data = self
            .col_iter()
            .flat_map(|col| col.iter().cloned().map(|c| c.norm()).collect::<Vec<_>>());
        let (nrows, ncols) = self.shape();
        Mat::from_fn(nrows, ncols, |_, _| col_wise_data.next().unwrap())
    }

    fn phase(&self) -> Option<Self::Output> {
        let mut col_wise_data = self
            .col_iter()
            .flat_map(|col| col.iter().map(|c| c.arg()).collect::<Vec<_>>());
        let (nrows, ncols) = self.shape();
        Some(Mat::from_fn(nrows, ncols, |_, _| {
            col_wise_data.next().unwrap()
        }))
    }
}
#[cfg(feature = "faer")]
impl Cartesian2Polar for Mat<f64> {
    type Output = Mat<f64>;

    fn magnitude(&self) -> Self::Output {
        self.clone()
    }

    fn phase(&self) -> Option<Self::Output> {
        None
    }
}

#[cfg(feature = "nalgebra")]
impl Cartesian2Polar for DMatrix<if64> {
    type Output = DMatrix<f64>;

    fn magnitude(&self) -> Self::Output {
        self.map(|x| x.modulus())
    }

    fn phase(&self) -> Option<Self::Output> {
        Some(self.map(|x| x.argument()))
    }
}

impl Cartesian2Polar for if64 {
    type Output = f64;

    fn magnitude(&self) -> Self::Output {
        self.norm()
    }

    fn phase(&self) -> Option<Self::Output> {
        Some(self.arg())
    }
}

pub trait Get {
    fn get(&self, row: usize, col: usize) -> &f64;
}
#[cfg(feature = "faer")]
impl Get for Mat<f64> {
    fn get(&self, row: usize, col: usize) -> &f64 {
        self.get(row, col)
    }
}
#[cfg(feature = "nalgebra")]
impl Get for DMatrix<f64> {
    fn get(&self, row: usize, col: usize) -> &f64 {
        &self[(row, col)]
    }
}
