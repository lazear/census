use std::u32;

/// Calaculate the sum of a slice
#[inline]
pub fn sum(slice: &[f64]) -> f64 {
    slice.iter().fold(0.0f64, |acc, x| acc + x)
}

/// Calculate the mean value of a slice
#[inline]
pub fn mean(slice: &[u32]) -> f64 {
    slice.iter().sum::<u32>() as f64 / slice.len() as f64
}

/// Return the maximum value of a slice
#[inline]
pub fn max(slice: &[u32]) -> u32 {
    slice
        .iter()
        .fold(u32::MIN, |acc, &x| if x > acc { x } else { acc })
}

/// Calculate the standard deviation (population) of a slice
#[inline]
pub fn stddev(slice: &[u32]) -> f64 {
    let mean = mean(slice);
    (slice
        .iter()
        .fold(0.0f64, |acc, x| acc + (*x as f64 - mean).powi(2))
        / slice.len() as f64)
        .sqrt()
}

/// Calculate the standard error (population) of a slice
#[inline]
pub fn stderr(slice: &[u32]) -> f64 {
    stddev(slice) / (slice.len() as f64).sqrt()
}

pub fn cv(slice: &[u32]) -> f64 {
    stddev(slice) / mean(slice)
}
