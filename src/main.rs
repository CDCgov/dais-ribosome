#![feature(let_chains)]
#![allow(unused_variables, dead_code)]

use clap::Parser;
//use rayon::{ThreadPoolBuilder, iter::ParallelBridge, prelude::ParallelIterator};
use std::path::PathBuf;
use zoe::{data::fasta::FastaSeq, prelude::*};

// If we later want to match the shell script, we can use:
// <https://docs.rs/git-version/latest/git_version/macro.git_describe.html>
const PROGRAM_VERSION: &str = concat!(env!("CARGO_PKG_NAME"), " v", env!("CARGO_PKG_VERSION"));

#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
/// CDS and and amino acid annotation tool for viruses.
pub struct Args {
    /// Data file to annotate in TSV or FASTA format.
    ///         †if classified, FASTA:  >ID|type_segment[_subtype]
    ///         *if classified, TSV:    ID<TAB>type_segment_[subtype]<TAB>sequence
    data_file: PathBuf,

    /// CDS and AA output, including coordinate mapping information, as a
    /// filename or path.
    sequence_output: Option<PathBuf>,

    /// Insertion output filename or path.
    insertion_output: Option<PathBuf>,

    /// Deletion output filename or path.
    deletion_output: Option<PathBuf>,

    /// Genomic file output prefix for sequences, insertion, and deletion.
    genomic_output_prefix: Option<PathBuf>,

    /// Name of the alignment module
    #[arg(short, long, default_value = "INFLUENZA")]
    module: String,

    /// Run in simultaneous multi-threaded mode.
    #[arg(short = 'T', long)]
    threads: bool,

    /// Write data as parquet files.
    #[arg(short = 'q', long)]
    output_parquet: bool,
}

fn main() {
    let args = Args::parse();

    println!("{mod}", mod=args.module);
}

struct Record {
    id: String,
    nucleotides: Nucleotides,
    classification: Option<String>,
}

impl TryFrom<FastaSeq> for Record {
    type Error = RibosomeError;
    fn try_from(datum: FastaSeq) -> Result<Self, Self::Error> {
        let FastaSeq { mut name, sequence } = datum;
        let classification = parse_id(&mut name);

        Ok(Record {
            id: name,
            nucleotides: sequence.into(),
            classification,
        })
    }
}

#[allow(dead_code)]
#[derive(Debug)]
enum RibosomeError {
    InvalidFastaFormat,
    InvalidTSVFormat,
    BlankFirstLIne,
}

impl std::fmt::Display for RibosomeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RibosomeError::InvalidFastaFormat => write!(f, "Invalid FASTA format"),
            RibosomeError::InvalidTSVFormat => write!(f, "Invalid TSV format"),
            RibosomeError::BlankFirstLIne => write!(f, "Blank first line"),
        }
    }
}

impl std::error::Error for RibosomeError {}

fn parse_id(id: &mut String) -> Option<String> {
    if let Some(offset) = id.find('|')
        && !id[offset + 1..].is_empty()
        && valid_classification(&id[offset + 1..])
    {
        let classification = id[offset + 1..].to_string();
        id.truncate(offset);
        Some(classification)
    } else {
        None
    }
}

// Rather than Regex, we will check for its existence in our module set
fn valid_classification(c: &str) -> bool {
    true
}

#[cfg(test)]
mod test {
    use crate::parse_id;

    #[test]
    fn parse_id_test() {
        let mut s = "ID|ANNOT".to_string();
        let annot = parse_id(&mut s);
        assert_eq!(("ID", Some("ANNOT".to_string())), (s.as_str(), annot));
    }
}
