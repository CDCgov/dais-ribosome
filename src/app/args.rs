use crate::app::io::{open_partition_writer, open_writer, optional_writers};
use clap::Parser;
use dais_ribosome::{error::RibosomeError, tsv::Writers};
use std::{fs::File, io::BufWriter, num::NonZero, path::PathBuf};

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
    #[arg(group = "output", requires_all = ["insertion_output", "deletion_output"])]
    sequence_output: Option<PathBuf>,

    /// Insertion output filename or path.
    #[arg(requires_all = ["sequence_output", "deletion_output"])]
    insertion_output: Option<PathBuf>,

    /// Deletion output filename or path.
    #[arg(requires_all = ["sequence_output", "insertion_output"])]
    deletion_output: Option<PathBuf>,

    /// Genomic file output prefix for sequences, insertion, and deletion.
    genomic_output_prefix: Option<PathBuf>,

    /// Name of the alignment module
    #[arg(short, long, default_value = "flu")]
    pub module: String,

    /// Run in simultaneous multi-threaded mode.
    #[arg(short = 'T', long)]
    pub threads: Option<NonZero<usize>>,

    // TODO: This does nothing
    /// Write data as parquet files.
    #[arg(short = 'q', long)]
    output_parquet: bool,

    /// Automatically detect the array task id from SGE or Slurm environment
    /// variables and write partition files for downstream collation.
    ///
    /// Output files are required and will be suffixed with a partition id.
    #[arg(short = 'G', long, requires = "output", conflicts_with_all = ["threads", "submit_grid_job"])]
    pub is_grid_task: bool,

    /// Submit and block on a grid engine (SGE or Slurm) array job of the
    /// specified size.
    #[arg(short = 'S', long, conflicts_with_all = ["threads", "is_grid_task"])]
    pub submit_grid_job: Option<usize>,

    /// Prints warning messages to stderr if any input sequences failed to be
    /// processed.
    #[arg(long, conflicts_with_all = ["is_grid_task", "submit_grid_job"])]
    pub verbose: bool,
}

impl Args {
    pub fn get_writers(&self) -> Result<Writers<BufWriter<File>>, RibosomeError> {
        Ok(Writers {
            seq: open_writer(&self.sequence_output, "seq")?,
            ins: open_writer(&self.insertion_output, "ins")?,
            del: open_writer(&self.deletion_output, "del")?,
        })
    }

    pub fn get_optional_writers(&self) -> Result<Option<Writers<BufWriter<File>>>, RibosomeError> {
        optional_writers(&self.genomic_output_prefix)
    }

    /// Return the output paths paired with their extension labels for grid
    /// partition collation.
    pub fn output_paths_for_grid(&self) -> Vec<(PathBuf, &str)> {
        let mut paths = Vec::new();
        if let Some(ref p) = self.sequence_output {
            paths.push((p.clone(), "seq"));
        }
        if let Some(ref p) = self.insertion_output {
            paths.push((p.clone(), "ins"));
        }
        if let Some(ref p) = self.deletion_output {
            paths.push((p.clone(), "del"));
        }
        // TODO: genomic output grid support
        paths
    }

    /// Build partition-suffixed writers for a grid array task.
    pub fn get_grid_writers(&self, task_id: usize) -> Result<Writers<BufWriter<File>>, RibosomeError> {
        Ok(Writers {
            seq: open_partition_writer(&self.sequence_output, task_id, "seq")?,
            ins: open_partition_writer(&self.insertion_output, task_id, "ins")?,
            del: open_partition_writer(&self.deletion_output, task_id, "del")?,
        })
    }
}
