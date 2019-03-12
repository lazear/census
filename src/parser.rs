//! Parse Census proteomics file
//!
//! Citation for census algorithm:
//!
//! Census 2: Isobaric labeling data analysis, Bioinformatics. 2014 Mar 28
//! A quantitative analysis software tool for mass spectrometry.based proteomics. Sung Kyu Park, John D Venable, Tao Xu, John R Yates III, Nature Methods, 2008, 5, 319-322

use super::*;

use std::iter::Peekable;
use std::str::Lines;

#[derive(PartialEq, PartialOrd, Debug)]
pub enum ErrorKind {
    /// Invalid beginning of line
    Invalid(char),
    /// Error converting to number
    Conversion,
    /// Unexpected end-of-file
    EOF,
}

/// Error that may occur during parsing of a Census file
#[derive(PartialEq, PartialOrd, Debug)]
pub struct Error {
    kind: ErrorKind,
    line: usize,
}

pub struct Parser<'s> {
    iter: Peekable<Lines<'s>>,
    /// Number of TMT channels to parse
    channels: u8,
    line: usize,
}

impl<'s> Parser<'s> {
    /// Create a new parser operating on input data
    pub fn new(input: &'s str) -> Parser<'s> {
        Parser {
            iter: input.lines().peekable(),
            channels: 0,
            line: 1,
        }
    }

    /// Convenience function for creating Error struct
    fn err(&self, kind: ErrorKind) -> Error {
        Error {
            kind,
            line: self.line,
        }
    }

    fn peek(&mut self) -> Option<&&'s str> {
        self.iter.peek()
    }

    fn next(&mut self) -> Option<&'s str> {
        let n = self.iter.next();
        if n.is_some() {
            self.line += 1;
        }
        n
    }

    fn parse_peptide(&mut self) -> Result<Peptide<'s>, Error> {
        let line = self.iter.next().ok_or_else(|| self.err(ErrorKind::EOF))?;
        let mut data = line.split_whitespace();
        assert_eq!(data.next(), Some("S"));
        let unique: bool = data.next().ok_or_else(|| self.err(ErrorKind::EOF))? == "U";
        let sequence = data.next().ok_or_else(|| self.err(ErrorKind::EOF))?;

        let mut values = Vec::with_capacity(self.channels as usize);

        for _ in 0..self.channels {
            let mz = data
                .next()
                .ok_or_else(|| self.err(ErrorKind::EOF))?
                .parse::<f64>()
                .map_err(|_| self.err(ErrorKind::Conversion))?;
            // discard normalized data
            let _ = data.next().ok_or_else(|| self.err(ErrorKind::EOF))?;
            values.push(mz);
        }
        Ok(Peptide {
            sequence,
            unique,
            values,
        })
    }

    fn parse_protein(&mut self) -> Result<Protein<'s>, Error> {
        let line = self.iter.next().ok_or_else(|| self.err(ErrorKind::EOF))?;
        let mut data = line.split('\t');
        assert_eq!(data.next(), Some("P"));
        let accession = data.next().ok_or_else(|| self.err(ErrorKind::EOF))?;
        let spectral_count = data
            .next()
            .ok_or_else(|| self.err(ErrorKind::EOF))?
            .parse::<u16>()
            .map_err(|_| self.err(ErrorKind::Conversion))?;
        let sequence_count = data
            .next()
            .ok_or_else(|| self.err(ErrorKind::EOF))?
            .parse::<u16>()
            .map_err(|_| self.err(ErrorKind::Conversion))?;
        let sequence_coverage = data
            .next()
            .ok_or_else(|| self.err(ErrorKind::EOF))?
            .trim_end_matches('%')
            .parse::<f32>()
            .map_err(|_| self.err(ErrorKind::Conversion))?;
        let molecular_weight = data
            .next()
            .ok_or_else(|| self.err(ErrorKind::EOF))?
            .parse::<u32>()
            .map_err(|_| self.err(ErrorKind::Conversion))?;

        let mut description = "";
        for n in data {
            description = n;
        }

        let mut peptides = Vec::new();
        while let Some(next) = self.iter.peek() {
            if next.starts_with('S') {
                peptides.push(self.parse_peptide()?);
            } else {
                // Next line should be a protein entry
                break;
            }
        }

        Ok(Protein {
            accession,
            spectral_count,
            sequence_count,
            sequence_coverage,
            molecular_weight,
            description,
            peptides,
        })
    }

    fn parse_headers(&mut self) -> Option<()> {
        while let Some(line) = self.peek() {
            if line.starts_with('H') {
                let line = self.next()?;
                if line.contains("m/z") {
                    self.channels = (line.matches("m/z_").count() / 2) as u8;
                }
            } else {
                return Some(());
            }
        }
        None
    }

    pub fn parse(mut self) -> Result<Dataset<'s>, Error> {
        let mut data = Vec::new();

        while let Some(line) = self.peek() {
            let init = line
                .chars()
                .next()
                .ok_or_else(|| self.err(ErrorKind::EOF))?;
            match init {
                'H' => self
                    .parse_headers()
                    .ok_or_else(|| self.err(ErrorKind::EOF))?,
                'P' => data.push(self.parse_protein()?),
                _ => return Err(self.err(ErrorKind::Invalid(init))),
            }
        }

        Ok(Dataset {
            proteins: data,
            channels: self.channels,
        })
    }
}
