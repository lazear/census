//! Parse Census proteomics file
//!

use super::*;

use std::fmt;
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

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Error parsing file at line {}: {:?}",
            self.line, self.kind
        )
    }
}

impl std::error::Error for Error {}

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

    fn parse_peptide(&mut self) -> Result<Peptide, Error> {
        let line = self.iter.next().ok_or_else(|| self.err(ErrorKind::EOF))?;
        // Using split_whitespace obfuscates missing 'U' values, and messes up
        // parsing
        let mut data = line.split('\t');
        assert_eq!(data.next(), Some("S"));

        let n = data.next().ok_or_else(|| self.err(ErrorKind::EOF))?;
        assert!(n.len() <= 1);
        let unique: bool = n == "U";
        let sequence = data.next().ok_or_else(|| self.err(ErrorKind::EOF))?.into();

        let mut values = Vec::with_capacity(self.channels as usize);

        for _ in 0..self.channels {
            let mz = data
                .next()
                .ok_or_else(|| self.err(ErrorKind::EOF))?
                .parse::<u32>()
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

    fn parse_protein(&mut self) -> Result<Protein, Error> {
        let line = self.iter.next().ok_or_else(|| self.err(ErrorKind::EOF))?;
        let mut data = line.split('\t');
        assert_eq!(data.next(), Some("P"));
        let accession = data.next().ok_or_else(|| self.err(ErrorKind::EOF))?.into();
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

        // let mut description = String::new();
        // for n in data {
        //     description = n.into();
        // }
        let description = data.last().ok_or_else(|| self.err(ErrorKind::EOF))?.into();

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
            channels: self.channels,
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

    pub fn parse(mut self) -> Result<Dataset, Error> {
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
