//! Utilities for filtering datasets and proteins based on a set of composable
//! rules
//!
//! ```ignore
//! protein:
//!     spectral_count: 10
//! peptide:
//!     exclude "C"
//!     tryptic
//!     intensity 5000
//!     unique
//!     channelcv 1, 6 10000.0
//!
//! ```
//!
use super::*;
use std::collections::HashSet;
mod parse;

/// Protein-level filter
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub enum ProteinFilter {
    /// Pass through proteins that have spectral counts >= N
    SpectralCounts(u16),
    /// Pass through proteins that have sequence counts >= N
    SequenceCounts(u16),
}

/// Peptide-level filter
///
/// Filter individual peptides within a protein based on sequence mactches,
/// total intensities, coefficient of variance between channels, or intensity
/// values on specified channels.
///
/// Peptide can also be filtered based on whether they have 2 tryptic ends,
/// or if they are unique.
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum PeptideFilter<'a> {
    /// Include only peptides that have a sequence matching the pattern
    SequenceMatch(&'a str),
    /// Exclude all peptides that have a sequence matching the pattern
    SequenceExclude(&'a str),
    /// Pass through peptides that have a total ion itensity >= N
    TotalIntensity(f64),

    /// ChannelCV(channels, N)
    ///
    /// Pass peptide if the coeff. of variance is < N between
    /// the specified channels
    ChannelCV(Vec<usize>, f64),

    /// ChannelIntensity(channel, cutoff)
    ///
    /// Pass through peptides that have an ion intensity >= N
    /// in the specified channel
    ChannelIntensity(usize, f64),

    /// Pass through tryptic peptides
    Tryptic,
    /// Include only unique peptides
    Unique,
}

/// Provides filtering functionality on datasets and proteins
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Filter<'a> {
    peptide_filters: Vec<PeptideFilter<'a>>,
    protein_filters: Vec<ProteinFilter>,
}

impl<'a> Default for Filter<'a> {
    fn default() -> Self {
        Filter {
            peptide_filters: Vec::new(),
            protein_filters: Vec::new(),
        }
    }
}

impl<'a> Filter<'a> {
    pub fn add_protein_filter(mut self, filter: ProteinFilter) -> Self {
        self.protein_filters.push(filter);
        self
    }

    pub fn add_peptide_filter(mut self, filter: PeptideFilter<'a>) -> Self {
        self.peptide_filters.push(filter);
        self
    }

    pub fn filter_dataset<'s>(&self, dataset: Dataset<'s>) -> Dataset<'s> {
        Dataset {
            channels: dataset.channels,
            proteins: dataset
                .proteins
                .into_iter()
                .filter_map(|prot| self.filter_protein(prot))
                .collect(),
        }
    }

    pub fn filter_protein<'s>(&self, mut protein: Protein<'s>) -> Option<Protein<'s>> {
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
                        if peptide.values.iter().fold(0.0, |acc, &x| acc + x) < *n {
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
                }
            }

            if pass {
                filtered.push(peptide)
            }
        }

        protein.peptides = filtered;

        let spec = protein.peptides.len() as u16;

        if spec == 0 {
            return None;
        }

        let seq = protein
            .peptides
            .iter()
            .map(|pep| pep.sequence)
            .collect::<HashSet<_>>()
            .len() as u16;

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
            }
        }

        Some(protein)
    }
}
