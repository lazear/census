use census::*;
use std::fs;

fn main() {
    let file = fs::read_to_string("./examples/data.txt").expect("Error opening file!");
    let data = census::read_census(&file).expect("Error parsing census file!");

    println!(
        "dataset containing {} proteins with {} channels",
        data.proteins.len(),
        data.channels
    );

    let filter = Filter::default()
        .add_protein_filter(ProteinFilter::SpectralCounts(2))
        .add_peptide_filter(PeptideFilter::ChannelIntensity(1, 5000.0))
        .add_peptide_filter(PeptideFilter::Unique)
        .add_peptide_filter(PeptideFilter::Tryptic)
        .add_peptide_filter(PeptideFilter::TotalIntensity(10_000.0));

    let data = data.filter(&filter);

    println!(
        "dataset containing {} proteins with {} channels",
        data.proteins.len(),
        data.channels
    );

    let first = data.proteins.first().unwrap();

    println!("{}, {}", first.accession, first.description);

    for peptide in &first.peptides {
        println!("{} {:?}", peptide.sequence, peptide.values);
    }
}
