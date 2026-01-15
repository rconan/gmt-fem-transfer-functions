use std::{fmt::Display, ops::Deref};

use serde::Serialize;

use super::{Cartesian2Polar, Dims, Get};

/// Frequency response data point
///
/// Frequency response magnitude and phase matrices at one frequency
#[derive(Debug, Serialize)]
pub struct FrequencyResponseData<T: Cartesian2Polar> {
    pub frequency: f64,
    pub magnitude: <T as Cartesian2Polar>::Output,
    pub phase: Option<<T as Cartesian2Polar>::Output>,
    pub u: Option<T>,
    pub v: Option<T>,
}
impl<T: Cartesian2Polar> FrequencyResponseData<T> {
    /// Creates a [FrequencyResponseData] instance from a frequency and response complex matrix
    pub fn new(frequency: f64, response: T) -> Self {
        Self {
            frequency,
            magnitude: response.magnitude(),
            phase: response.phase(),
            u: None,
            v: None,
        }
    }
    pub fn new_svd(frequency: f64, (response, u, v): (T, Option<T>, Option<T>)) -> Self {
        Self {
            frequency,
            magnitude: response.magnitude(),
            phase: response.phase(),
            u,
            v,
        }
    }
}
impl<T> Display for FrequencyResponseData<T>
where
    T: Cartesian2Polar,
    <T as Cartesian2Polar>::Output: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(phase) = &self.phase {
            write!(f, "{},{},{}", self.frequency, self.magnitude, phase)
        } else {
            write!(f, "{},{}", self.frequency, self.magnitude)
        }
    }
}

/// Collection of [FrequencyResponseData]
#[derive(Debug, Serialize)]
pub struct FrequencyResponseVec<T: Cartesian2Polar>(
    #[serde(rename = "data")] Vec<FrequencyResponseData<T>>,
);
impl<T: Cartesian2Polar> Default for FrequencyResponseVec<T> {
    fn default() -> Self {
        Self(vec![])
    }
}

/// Frequency response extremum
pub struct Extremum {
    pub i: usize,
    pub x: f64,
    pub y: f64,
}

impl<T> FrequencyResponseVec<T>
where
    T: Cartesian2Polar,
    <T as Cartesian2Polar>::Output: Get,
{
    /// Creates a new [FrequencyResponseVec] instance from a vector of [FrequencyResponseData]
    pub fn new(frequency_response_datas: Vec<FrequencyResponseData<T>>) -> Self {
        Self(frequency_response_datas)
    }
    /// Returns the frequency vector
    pub fn frequencies(&self) -> Vec<f64> {
        self.iter().map(|fr| fr.frequency).collect()
    }
    /// Returns the extrema of the frequency response
    pub fn extrema(&self, indices: Option<(usize, usize)>) -> Vec<Extremum> {
        let (row, col) = indices.unwrap_or_default();
        let mut iter = self
            .0
            .iter()
            .map(|data| (data.frequency, data.magnitude.get(row, col)))
            .peekable();
        let mut previous_s = Option::<f64>::None;
        let mut extrema = vec![];
        let mut i = 0;
        while let Some(((fa, ma), &(_, mb))) = iter.next().zip(iter.peek()) {
            let s = (mb - ma).signum();
            if let Some(previous_s) = previous_s
                && previous_s != s
            {
                extrema.push(Extremum { i, x: fa, y: *ma })
            }
            previous_s = Some(s);
            i += 1;
        }
        extrema
    }
}

impl<T: Cartesian2Polar> Deref for FrequencyResponseVec<T> {
    type Target = [FrequencyResponseData<T>];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Cartesian2Polar> FromIterator<FrequencyResponseData<T>> for FrequencyResponseVec<T> {
    fn from_iter<I: IntoIterator<Item = FrequencyResponseData<T>>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<T> Display for FrequencyResponseVec<T>
where
    T: Cartesian2Polar,
    <T as Cartesian2Polar>::Output: Get,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "GMT FEM frequency response matrix {}x{:?}",
            self.len(),
            self[0].magnitude.size()
        )?;
        match self.len() {
            n if n == 1 => {
                writeln!(f, " @ {:.2}Hz", self[0].frequency)
            }
            n if n < 6 => {
                writeln!(f, " @ {:.2?}Hz", self.frequencies())
            }
            _ => {
                writeln!(
                    f,
                    " @ [{:.2},{:.2}]Hz",
                    self[0].frequency,
                    self.last().unwrap().frequency
                )
            }
        }
    }
}
