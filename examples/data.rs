//! This is an example of how to use the census crate to read in some sample
//! data (not supplied, sorry.), parse, filter, and normalize the data to
//! the ion intensities from channel 1. The normalized data is then written
//! to a new output file
//!
//! Also demonstrated here is the ability to use the "serialization" feature
//! which allows `Filter` objects to read or written from a file in JSON
//! format.  
//!
//! This whole process should execute in ~250ms for a 25Mb file of raw data.
use census_proteomics::*;
#[cfg(feature = "serialization")]
use serde_json;
use std::fs;
use std::io::prelude::*;

fn main() -> std::io::Result<()> {
    let file = fs::read_to_string("./examples/data.txt")?;
    let data = census_proteomics::read_census(&file).expect("Error parsing census file!");

    #[cfg(feature = "serialization")]
    let s = fs::read_to_string("./examples/filter.json")?;

    #[cfg(feature = "serialization")]
    let filter = match serde_json::from_str(&s) {
        Ok(f) => f,
        Err(e) => {
            println!("Error while parsing filter.json {:?}", e);
            std::process::abort();
        }
    };

    #[cfg(not(feature = "serialization"))]
    let filter = Filter::default()
        .add_peptide_filter(PeptideFilter::ChannelIntensity(1, 1000))
        .add_peptide_filter(PeptideFilter::Unique)
        .add_peptide_filter(PeptideFilter::Tryptic)
        .add_peptide_filter(PeptideFilter::TotalIntensity(5000))
        .add_peptide_filter(PeptideFilter::Purity(0.9));

    let data = data.filter(&filter);
    let mut output = fs::File::create("out.txt")?;

    writeln!(
        output,
        "accession\tdescription\tspectral_count\tsequence_count\tsequence\tscan\tpurity\t{}",
        (1..=data.channels)
            .map(|i| format!("channel_{}", i))
            .collect::<Vec<String>>()
            .join("\t")
    )?;

    for prot in &data.proteins {
        for pep in &prot.peptides {
            // Write our data to a tab-delimited file
            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                prot.accession,
                prot.description,
                prot.spectral_count,
                prot.sequence_count,
                pep.sequence,
                pep.scan,
                pep.purity,
                pep.values
                    .iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<String>>()
                    .join("\t")
            )?;
        }
    }

    Ok(())
}
