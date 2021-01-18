#[cfg(feature = "serialization")]
use serde::Serialize;

#[cfg_attr(feature = "serde", derive(Serialize))]
#[derive(PartialEq, PartialOrd, Clone, Default)]
/// Protein-level TMT quantification data, as well as additional
/// metadata about the protein that is output in the Census file
pub struct Protein {
    /// Uniprot accession identifier
    pub accession: String,
    /// Long-form description
    pub description: String,
    /// Number of spectral counts
    pub spectral_count: u16,
    /// Number of unique sequence counts
    pub sequence_count: u16,
    /// Sequence coverage
    pub sequence_coverage: f32,
    /// Molecular weight
    pub molecular_weight: u32,
    /// Raw signal intensity channels
    pub peptides: Vec<Peptide>,

    pub channels: u8,
}

impl Protein {
    /// Return the summed intensities for all peptides
    pub fn total(&self) -> Vec<u32> {
        let mut v = Vec::with_capacity(self.channels as usize);
        for c in 0..self.channels {
            let sum = self.peptides.iter().map(|pep| pep.values[c as usize]).sum();
            v.push(sum);
        }
        v
    }

    /// Return a vector of normalized ratios, where the signal intensity
    /// for each channel is divided by the sum of all channels
    pub fn ratios(&self) -> Vec<f64> {
        let values = self.total();
        let total = values.iter().sum::<u32>() as f64;
        values.iter().map(|v| *v as f64 / total).collect()
    }
}

#[cfg_attr(feature = "serde", derive(Serialize))]
#[derive(PartialEq, PartialOrd, Clone, Debug, Default)]
/// Peptide-level TMT quantification data
pub struct Peptide {
    /// Peptide sequence
    pub sequence: String,
    /// Raw isobaric ion intensity values
    pub values: Vec<u32>,
    /// Is this a unique peptide?
    pub unique: bool,

    pub purity: f32,

    pub scan: usize,
}

impl Peptide {
    /// Return a boolean indicating whether the peptide has 2 tryptic sites
    pub fn tryptic(&self) -> bool {
        let cterm = self.sequence.ends_with('-');
        let front = self.sequence.starts_with(|c| match c {
            'K' | 'R' | '-' => true,
            _ => false,
        });
        let end = self
            .sequence
            .split('.')
            .skip(1)
            .next()
            .map(|s| {
                s.ends_with(|c| match c {
                    'K' | 'R' => true,
                    _ => cterm,
                })
            })
            .unwrap_or(false);
        front && end
        // let kdot = self.sequence.matches("K.").count();
        // let rdot = self.sequence.matches("R.").count();
        // kdot + rdot == 2
    }

    /// Return a vector of normalized ratios, where the signal intensity
    /// for each channel is divided by the sum of all channels
    pub fn ratios(&self) -> Vec<f64> {
        let total: f64 = self.values.iter().sum::<u32>() as f64;
        self.values.iter().map(|v| *v as f64 / total).collect()
    }

    /// Swap channels A and B, which are 0 indexed into the peptide values
    /// vector.
    ///
    /// # May panic
    ///
    /// May panic if A or B exceed the length of the vector
    pub fn swap_channels(&mut self, a: usize, b: usize) {
        self.values.swap(a, b)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn gen_peptide(sequence: &str) -> Peptide {
        Peptide {
            sequence: sequence.into(),
            ..Peptide::default()
        }
    }
    #[test]
    fn test_trypic() {
        assert!(gen_peptide("-.KMDKDK.-").tryptic());
        assert!(!gen_peptide("S.KMDKDK.-").tryptic());
        assert!(gen_peptide("R.KMDKDK.-").tryptic());
        assert!(!gen_peptide("K.KMDKDT.A").tryptic());
    }
}
