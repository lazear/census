//! Utilities for filtering datasets and proteins based on a set of composable
//! rules
use super::*;
#[cfg(feature = "serialization")]
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

/// Protein-level filter
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub enum ProteinFilter {
    /// Include only proteins that have spectral counts >= N
    SpectralCounts(u16),
    /// Include only proteins that have sequence counts >= N
    SequenceCounts(u16),
    /// Include only proteins that do not have "Reverse" in their
    /// UniProt accession
    ExcludeReverse,
}

/// Peptide-level filter
///
/// Filter individual peptides within a protein based on sequence mactches,
/// total intensities, coefficient of variance between channels, or intensity
/// values on specified channels.
///
/// Peptide can also be filtered based on whether they have 2 tryptic ends,
/// or if they are unique.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum PeptideFilter<'a> {
    /// Include only peptides that have a sequence matching the pattern
    SequenceMatch(&'a str),
    /// Include only peptides that do NOT have a sequence matching the pattern
    SequenceExclude(&'a str),
    /// Include only peptides that have a total ion itensity >= N
    TotalIntensity(u32),

    /// Include only peptides where the total intensity in a set of channels
    /// >= N
    TotalIntensityChannels(Vec<usize>, u32),

    /// ChannelCV(channels, N)
    ///
    /// Include only peptides where the coeff. of variance is < N between
    /// the specified channels
    ChannelCV(Vec<usize>, f64),

    /// ChannelIntensity(channel, cutoff)
    ///
    /// Include only peptides that have an ion intensity >= N
    /// in the specified channel
    ChannelIntensity(usize, u32),

    /// Include only tryptic peptides
    Tryptic,
    /// Include only unique peptides
    Unique,
}

/// Provides filtering functionality on datasets and proteins
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Filter<'a> {
    #[cfg_attr(feature = "serde", serde(borrow))]
    peptide_filters: Vec<PeptideFilter<'a>>,
    protein_filters: Vec<ProteinFilter>,
}

impl<'a> Default for Filter<'a> {
    /// Construct a filter with no rules
    fn default() -> Self {
        Filter {
            peptide_filters: Vec::new(),
            protein_filters: Vec::new(),
        }
    }
}

impl<'a> Filter<'a> {
    /// Add a new `ProteinFilter` to the `Filter` object.
    ///
    /// This follows the Builder pattern
    pub fn add_protein_filter(mut self, filter: ProteinFilter) -> Self {
        self.protein_filters.push(filter);
        self
    }

    /// Add a new `PeptideFilter` to the `Filter` object.
    ///
    /// This follows the Builder pattern
    pub fn add_peptide_filter(mut self, filter: PeptideFilter<'a>) -> Self {
        self.peptide_filters.push(filter);
        self
    }

    /// Return a new `Dataset` that only contains filtered `Protein`'s
    pub fn filter_dataset(&self, dataset: Dataset) -> Dataset {
        Dataset {
            channels: dataset.channels,
            proteins: dataset
                .proteins
                .into_iter()
                .filter_map(|prot| self.filter_protein(prot))
                .collect(),
        }
    }

    /// Filter a `Protein`, returning `Some` if it passes any
    /// `ProteinFilter`s that need to be applied or `None` if the protein
    /// fails a given `ProteinFilter`.
    ///
    /// The peptides associated with the returned `Protein` object aree those
    /// that passed any given `PeptideFilter`s.
    pub fn filter_protein(&self, mut protein: Protein) -> Option<Protein> {
        // First run through any protein level filters
        for filter in &self.protein_filters {
            match filter {
                ProteinFilter::SequenceCounts(n) => {
                    if protein.sequence_count < *n {
                        return None;
                    }
                }
                ProteinFilter::SpectralCounts(n) => {
                    if protein.spectral_count < *n {
                        return None;
                    }
                }
                ProteinFilter::ExcludeReverse => {
                    if protein.accession.contains("Reverse") {
                        return None;
                    }
                }
            }
        }

        let mut filtered = Vec::new();

        // Iterate through all of the peptides in the protein container,
        // applying relevant filters as we go.
        for peptide in protein.peptides {
            let mut pass = true;
            for filter in &self.peptide_filters {
                match filter {
                    PeptideFilter::SequenceExclude(pat) => {
                        if peptide.sequence.contains(pat) {
                            pass = false;
                        }
                    }
                    PeptideFilter::SequenceMatch(pat) => {
                        if !peptide.sequence.contains(pat) {
                            pass = false;
                        }
                    }
                    PeptideFilter::TotalIntensity(n) => {
                        if peptide.values.iter().sum::<u32>() < *n {
                            pass = false;
                        }
                    }
                    PeptideFilter::Tryptic => {
                        if !peptide.tryptic() {
                            pass = false;
                        }
                    }
                    PeptideFilter::Unique => {
                        if !peptide.unique {
                            pass = false;
                        }
                    }
                    PeptideFilter::ChannelCV(channels, cutoff) => {
                        let mut v = Vec::new();
                        for chan in channels.iter() {
                            if chan - 1 < peptide.values.len() {
                                v.push(peptide.values[chan - 1]);
                            }
                        }
                        if util::cv(&v) >= *cutoff {
                            pass = false;
                        }
                    }
                    PeptideFilter::ChannelIntensity(channel, cutoff) => {
                        // Ignore incorrect channel values
                        if channel - 1 < peptide.values.len()
                            && peptide.values[channel - 1] < *cutoff
                        {
                            pass = false;
                        }
                    }
                    PeptideFilter::TotalIntensityChannels(chan, cutoff) => {
                        let mut sum = 0;
                        for c in chan {
                            if c - 1 < peptide.values.len() {
                                sum += peptide.values[*c - 1];
                            }
                        }
                        if sum < *cutoff {
                            pass = false;
                        }
                    }
                }
            }

            if pass {
                filtered.push(peptide)
            }
        }
        // We must have at least a single peptide...
        if filtered.len() == 0 {
            return None;
        }

        protein.peptides = filtered;

        let spec = protein.peptides.len() as u16;
        let seq = protein
            .peptides
            .iter()
            .map(|pep| &pep.sequence)
            .collect::<HashSet<_>>()
            .len() as u16;

        protein.spectral_count = spec;
        protein.sequence_count = seq;

        // Second pass through protein filters, in case we no longer have
        // enough filtered peptides
        for filter in &self.protein_filters {
            match filter {
                ProteinFilter::SequenceCounts(n) => {
                    if seq < *n {
                        return None;
                    }
                }
                ProteinFilter::SpectralCounts(n) => {
                    if spec < *n {
                        return None;
                    }
                }
                _ => {}
            }
        }

        Some(protein)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn total_intensity_channels() {
        let p1 = Peptide {
            sequence: "aa".into(),
            values: vec![1, 2998, 5000, 84, 4738, 9384],
            unique: true,
            scan: 0,
        };
        let p2 = Peptide {
            sequence: "aaa".into(),
            values: vec![10000, 0, 433, 61346, 41, 5555],
            unique: true,
            scan: 0,
        };

        let p3 = Peptide {
            sequence: "aaaa".into(),
            values: vec![1, 2999, 0, 0, 0, 0],
            unique: true,
            scan: 0,
        };

        let prot = Protein {
            accession: "".into(),
            description: "".into(),
            spectral_count: 10,
            sequence_count: 3,
            sequence_coverage: 0.3,
            molecular_weight: 10,
            peptides: vec![p1.clone(), p2.clone(), p3.clone()],
            channels: 6,
        };

        let mut fil = Filter {
            peptide_filters: Vec::new(),
            protein_filters: Vec::new(),
        };

        fil = fil.add_peptide_filter(PeptideFilter::TotalIntensityChannels(vec![1, 2], 3000));
        let p = fil.filter_protein(prot).unwrap();
        assert_eq!(p.peptides.len(), 2);
        assert_eq!(p.sequence_count, 2);
        assert_eq!(p.peptides, vec![p2.clone(), p3.clone()]);
    }
}
