use std::f64;

/// Calaculate the sum of a slice
#[inline]
pub fn sum(slice: &[f64]) -> f64 {
    slice.iter().fold(0.0f64, |acc, x| acc + x)
}

/// Calculate the mean value of a slice
#[inline]
pub fn mean(slice: &[f64]) -> f64 {
    sum(slice) / slice.len() as f64
}

/// Return the maximum value of a slice
#[inline]
pub fn max(slice: &[f64]) -> f64 {
    slice
        .iter()
        .fold(f64::MIN, |acc, &x| if x > acc { x } else { acc })
}

/// Calculate the standard deviation (population) of a slice
#[inline]
pub fn stddev(slice: &[f64]) -> f64 {
    let mean = mean(slice);
    (slice.iter().fold(0.0f64, |acc, x| acc + (x - mean).powi(2)) / slice.len() as f64).sqrt()
}

/// Calculate the standard error (population) of a slice
#[inline]
pub fn stderr(slice: &[f64]) -> f64 {
    stddev(slice) / (slice.len() as f64).sqrt()
}

pub fn cv(slice: &[f64]) -> f64 {
    stddev(slice) / mean(slice)
}
