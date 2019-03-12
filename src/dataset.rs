//! Collection of `Protein` objects representing a single dataset
use super::*;
use std::collections::{HashMap, HashSet};

/// Container for proteomics data read from a Census version file
pub struct Dataset<'s> {
    /// TMT data for each protein in the dataset
    pub proteins: Vec<Protein<'s>>,
    /// Number of TMT channels in the dataset
    pub channels: u8,
}

impl<'s> Dataset<'s> {
    /// Return a set of all UniProt KB accession ID's present in the
    /// `Dataset`
    pub fn accessions(&self) -> HashSet<&'s str> {
        self.proteins.iter().map(|pr| pr.accession).collect()
    }

    /// Create a `HashMap` correlating a UniProt KB accession to Protein-level
    /// quant data
    pub fn map(&self) -> HashMap<&'s str, &Protein<'s>> {
        self.proteins.iter().map(|pr| (pr.accession, pr)).collect()
    }

    pub fn filter(self, filter: &Filter) -> Self {
        filter.filter_dataset(self)
    }
}
