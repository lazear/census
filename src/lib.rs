//! A high-performance Rust library for parsing, filtering, and manipulating
//! multiplexed isobaric data that has been quantified using the Census
//! algorithm
mod dataset;
mod filter;
mod parser;
mod protein;
mod util;

pub use dataset::Dataset;
pub use filter::{Filter, PeptideFilter, ProteinFilter};
pub use parser::{Error, Parser};
pub use protein::{Peptide, Protein};

/// Parse a string containing a complete census file into a `Dataset`
pub fn read_census(input: &str) -> Result<Dataset, Error> {
    Parser::new(input).parse()
}
