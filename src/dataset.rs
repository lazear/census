//! Collection of `Protein` objects representing a single dataset
use super::*;
#[cfg(feature = "serialization")]
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

#[cfg_attr(feature = "serde", derive(Serialize))]
/// Container for proteomics data read from a Census version file
pub struct Dataset {
    /// TMT data for each protein in the dataset
    pub proteins: Vec<Protein>,
    /// Number of TMT channels in the dataset
    pub channels: u8,
}

impl Dataset {
    /// Return a set of all UniProt KB accession ID's present in the
    /// `Dataset`
    pub fn accessions(&self) -> HashSet<&'_ str> {
        self.proteins
            .iter()
            .map(|pr| pr.accession.as_ref())
            .collect()
    }

    /// Create a `HashMap` correlating a UniProt KB accession to Protein-level
    /// quant data
    pub fn map(&self) -> HashMap<&'_ str, &Protein> {
        self.proteins
            .iter()
            .map(|pr| (pr.accession.as_ref(), pr))
            .collect()
    }

    pub fn filter(self, filter: &Filter) -> Self {
        filter.filter_dataset(self)
    }
}
