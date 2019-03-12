//! Filter rule parser. Not very pretty code, but it's fast and it doesn't allocate
//! (except for in ChannelCV, where it allocates a Vec)
use super::*;

#[derive(PartialEq, PartialOrd, Debug)]
enum ErrorKind<'s> {
    Command(&'s str),
    Expected(&'s str),
    Conversion,
}

fn take_while<'s, F: Fn(char) -> bool>(input: &'s str, pred: F) -> (&'s str, &'s str) {
    let mut iter = input.chars();
    loop {
        let rest = iter.as_str();
        match iter.next() {
            Some(c) if pred(c) => {}
            _ => {
                if rest.len() != input.len() {
                    let ret = &input[..input.len() - rest.len()];
                    return (ret, rest);
                } else {
                    return ("", rest);
                }
            }
        }
    }
}

#[inline]
fn take_whitespace(input: &str) -> &str {
    input.trim_start_matches(char::is_whitespace)
}

#[inline]
fn take_word(input: &str) -> (&str, &str) {
    take_while(input, |ch| !ch.is_whitespace())
}

fn expect<'s>(input: &'s str, m: &'static str) -> Option<&'s str> {
    if input.starts_with(m) {
        Some(input.trim_start_matches(m))
    } else {
        None
    }
}

fn parse_protein_filter<'s>(mut input: &'s str) -> Result<(ProteinFilter, &'s str), ErrorKind> {
    let (cmd, input) = take_word(take_whitespace(input));
    match cmd.as_ref() {
        "spectral_counts" => {
            let input = take_whitespace(
                expect(take_whitespace(input), "=").ok_or(ErrorKind::Expected("="))?,
            );
            let (num, rest) = take_word(input);
            Ok((
                ProteinFilter::SpectralCounts(
                    num.parse::<u16>().map_err(|_| ErrorKind::Conversion)?,
                ),
                rest,
            ))
        }
        "sequence_counts" => {
            let input = take_whitespace(
                expect(take_whitespace(input), "=").ok_or(ErrorKind::Expected("="))?,
            );
            let (num, rest) = take_word(input);
            Ok((
                ProteinFilter::SequenceCounts(
                    num.parse::<u16>().map_err(|_| ErrorKind::Conversion)?,
                ),
                rest,
            ))
        }
        _ => Err(ErrorKind::Command(cmd)),
    }
}

fn parse_peptide_filter<'s>(mut input: &'s str) -> Result<(PeptideFilter, &'s str), ErrorKind> {
    let (cmd, input) = take_word(take_whitespace(input));
    match cmd.as_ref() {
        "total_intensity" => {
            let input = take_whitespace(
                expect(take_whitespace(input), "=").ok_or(ErrorKind::Expected("="))?,
            );
            let (num, rest) = take_word(input);
            Ok((
                PeptideFilter::TotalIntensity(
                    num.parse::<f64>().map_err(|_| ErrorKind::Conversion)?,
                ),
                rest,
            ))
        }
        "channel_intensity" => {
            let input = take_whitespace(
                expect(take_whitespace(input), "=").ok_or(ErrorKind::Expected("="))?,
            );
            let (chan, rest) = take_word(input);
            let (cutoff, rest) = take_word(take_whitespace(rest));

            Ok((
                PeptideFilter::ChannelIntensity(
                    chan.parse::<usize>().map_err(|_| ErrorKind::Conversion)?,
                    cutoff.parse::<f64>().map_err(|_| ErrorKind::Conversion)?,
                ),
                rest,
            ))
        }
        "channel_cv" => {
            // This is UGLY

            let mut channels = vec![];
            let mut input = take_whitespace(
                expect(take_whitespace(input), "=").ok_or(ErrorKind::Expected("="))?,
            );

            // Read the comma delimited channels
            loop {
                let (chan, rest) = take_while(take_whitespace(input), char::is_numeric);
                let rest = take_whitespace(rest);
                if !rest.trim_start().starts_with(",") {
                    break;
                }
                input = rest.trim_start_matches(",");
                channels.push(chan.parse::<usize>().map_err(|_| ErrorKind::Conversion)?);
            }

            // Read the last channel, it doesn't have a comma after it
            let (last, rest) = take_word(take_whitespace(input));

            channels.push(last.parse::<usize>().map_err(|_| ErrorKind::Conversion)?);

            // Now read the cutoff value
            let (cutoff, rest) = take_word(take_whitespace(rest));
            Ok((
                PeptideFilter::ChannelCV(
                    channels,
                    cutoff.parse::<f64>().map_err(|_| ErrorKind::Conversion)?,
                ),
                rest,
            ))
        }
        "sequence_match" => {
            let input = take_whitespace(
                expect(take_whitespace(input), "=").ok_or(ErrorKind::Expected("="))?,
            );
            let (seq, rest) = take_word(input);
            Ok((PeptideFilter::SequenceMatch(seq), rest))
        }
        "sequence_exclude" => {
            let input = take_whitespace(
                expect(take_whitespace(input), "=").ok_or(ErrorKind::Expected("="))?,
            );
            let (seq, rest) = take_word(input);
            Ok((PeptideFilter::SequenceExclude(seq), rest))
        }
        "tryptic" => Ok((PeptideFilter::Tryptic, input)),
        "unique" => Ok((PeptideFilter::Unique, input)),
        _ => Err(ErrorKind::Command(cmd)),
    }
}

fn parse<'s>(mut input: &'s str) -> Result<Filter, ErrorKind> {
    let mut filter = Filter::default();
    loop {
        dbg!(&input);
        if input.starts_with("protein") {
            let (_, mut sr) = take_word(input);
            loop {
                sr = sr.trim_start_matches(char::is_whitespace);
                if sr.starts_with("peptide") || sr.is_empty() {
                    input = sr;
                    break;
                }
                match parse_protein_filter(sr) {
                    Ok((filt, rest)) => {
                        filter = filter.add_protein_filter(filt);
                        sr = rest;
                    }
                    Err(e) => return Err(e),
                }
            }
        } else if input.starts_with("peptide") {
            let (_, mut sr) = take_word(input);
            loop {
                sr = sr.trim_start_matches(char::is_whitespace);
                if sr.starts_with("protein") || sr.is_empty() {
                    input = sr;
                    break;
                }
                match parse_peptide_filter(sr) {
                    Ok((filt, rest)) => {
                        filter = filter.add_peptide_filter(filt);
                        sr = rest;
                    }
                    Err(e) => return Err(e),
                }
            }
        } else {
            break;
        }
    }
    Ok(filter)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn protein_filters() {
        let input = "
        spectral_counts = 4
        sequence_counts = 2";
        let (fst, rest) = parse_protein_filter(input).unwrap();
        dbg!(&rest);
        let (snd, rest) = parse_protein_filter(rest).unwrap();
        assert_eq!(fst, ProteinFilter::SpectralCounts(4));
        assert_eq!(snd, ProteinFilter::SequenceCounts(2));

        assert!(rest.is_empty())
    }

    #[test]
    fn peptide_filters() {
        let input = "
        channel_cv = 1,2,3,4 0.5
        sequence_match = CCS tryptic    \t \n \r unique";
        let (fst, rest) = parse_peptide_filter(input).unwrap();
        let (snd, rest) = parse_peptide_filter(rest).unwrap();
        assert_eq!(fst, PeptideFilter::ChannelCV(vec![1, 2, 3, 4], 0.5));
        assert_eq!(snd, PeptideFilter::SequenceMatch("CCS"));
        let (thd, rest) = parse_peptide_filter(rest).unwrap();
        let (frt, rest) = parse_peptide_filter(rest).unwrap();
        assert_eq!(thd, PeptideFilter::Tryptic);
        assert_eq!(frt, PeptideFilter::Unique);

        assert!(rest.is_empty())
    }

    #[test]
    fn config() {
        let input = "protein:
    spectral_counts = 10
    sequence_counts = 2
peptide:
    channel_cv = 1, 2, 6,   7 0.05
    sequence_exclude = C
    tryptic
    unique";

        let expected = Filter::default()
            .add_peptide_filter(PeptideFilter::ChannelCV(vec![1, 2, 6, 7], 0.05))
            .add_peptide_filter(PeptideFilter::SequenceExclude("C"))
            .add_peptide_filter(PeptideFilter::Tryptic)
            .add_peptide_filter(PeptideFilter::Unique)
            .add_protein_filter(ProteinFilter::SpectralCounts(10))
            .add_protein_filter(ProteinFilter::SequenceCounts(2));

        assert_eq!(parse(input), Ok(expected));
    }
}
