use census::*;
#[cfg(feature = "serialization")]
use serde_json;
use std::fs;
use std::io::prelude::*;

fn main() -> std::io::Result<()> {
    let file = fs::read_to_string("./examples/data.txt")?;
    let data = census::read_census(&file).expect("Error parsing census file!");

    #[cfg(feature = "serialization")]
    let s = fs::read_to_string("./examples/filter.json")?;
    #[cfg(feature = "serialization")]
    let filter = serde_json::from_str(&s).unwrap();

    #[cfg(not(feature = "serialization"))]
    let filter = Filter::default()
        .add_peptide_filter(PeptideFilter::ChannelIntensity(1, 1000))
        .add_peptide_filter(PeptideFilter::Unique)
        .add_peptide_filter(PeptideFilter::Tryptic)
        .add_peptide_filter(PeptideFilter::TotalIntensity(5000));

    let data = data.filter(&filter);
    let mut output = fs::File::create("output.txt")?;

    writeln!(
        output,
        "accession\tdescription\tspectral_count\tsequence_count\t{}",
        (1..=data.channels)
            .map(|i| format!("channel_{}", i))
            .collect::<Vec<String>>()
            .join("\t")
    )?;

    for prot in &data.proteins {
        // Sum up all of the intensities for the constituent peptides
        // for this protein
        let vals = prot.total();
        // Take the average of the first 4 TMT channels
        let ctrl = vals.iter().take(4).sum::<u32>() as f64 / 4.0;
        // Divide each channel's intensity by the average intensity
        // of channels 1 - 4, and convert it to a tab delimited
        // string so we can write it.
        let adj = vals
            .into_iter()
            .map(|v| format!("{}", (v as f64) / ctrl))
            .collect::<Vec<String>>()
            .join("\t");

        // Write our data to a tab-delimited file
        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}",
            prot.accession, prot.description, prot.spectral_count, prot.sequence_count, adj
        )?;
    }

    Ok(())
}
