use super::io;
use clap::Parser;
use dais_ribosome::data::RibosomeError;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
/// CDS and and amino acid annotation tool for viruses.
pub struct Args {
    /// Data file to annotate in TSV or FASTA format.
    ///         `†if classified, FASTA:  >ID|type_segment[_subtype]`
    ///         `*if classified, TSV:    ID<TAB>type_segment_[subtype]<TAB>sequence`
    pub data_file: PathBuf,

    /// CDS and AA output, including coordinate mapping information, as a
    /// filename or path.
    #[arg(group = "output", requires_all = ["insertion_output","deletion_output"])]
    sequence_output: Option<PathBuf>,

    /// Insertion output filename or path.
    #[arg(requires_all = ["sequence_output","deletion_output"])]
    insertion_output: Option<PathBuf>,

    /// Deletion output filename or path.
    #[arg(requires_all = ["sequence_output","insertion_output"])]
    deletion_output: Option<PathBuf>,

    /// Genomic file output prefix for sequences, insertion, and deletion.
    genomic_output_prefix: Option<PathBuf>,

    /// Name of the alignment module
    #[arg(short, long, default_value = "flu")]
    pub module: String,

    /// Run in simultaneous multi-threaded mode.
    #[arg(short = 'T', long)]
    threads: bool,

    /// Write data as parquet files.
    #[arg(short = 'q', long)]
    output_parquet: bool,
}

impl Args {
    pub fn get_writers(&self) -> Result<io::Writers, RibosomeError> {
        Ok(io::Writers {
            seq: io::open_writer(&self.sequence_output, "seq")?,
            ins: io::open_writer(&self.insertion_output, "ins")?,
            del: io::open_writer(&self.deletion_output, "del")?,
        })
    }

    pub fn get_optional_writers(&self) -> Result<Option<io::Writers>, RibosomeError> {
        io::optional_writers(&self.genomic_output_prefix)
    }
}
