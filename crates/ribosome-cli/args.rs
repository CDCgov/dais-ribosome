use crate::{
    log::time_stamp,
    par_utils::grid::{GridCompatibleArgs, GridCompatibleCli},
};
use clap::{Arg, Parser, builder::OsStr};
use std::{
    num::NonZero,
    path::{Path, PathBuf},
};
use zoe::prelude::rand_sequence;

#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
#[command(after_help = "† For classified FASTA:  >ID|ctype\n  and classified TSV:    ID<TAB>ctype<TAB>sequence")]
/// CDS and and amino acid annotation tool for viruses.
pub struct Cli {
    /// Data file to annotate in TSV or FASTA format.†
    pub data_file: PathBuf,

    /// CDS and AA output, including coordinate mapping information, as a
    /// filename or path.
    #[arg(requires_all = ["insertion_output", "deletion_output"])]
    pub sequence_output: Option<PathBuf>,

    /// Insertion output filename or path.
    #[arg(requires_all = ["sequence_output", "deletion_output"])]
    pub insertion_output: Option<PathBuf>,

    /// Deletion output filename or path.
    #[arg(requires_all = ["sequence_output", "insertion_output"])]
    pub deletion_output: Option<PathBuf>,

    /// Genome sequence, insertion, and deletion output paths.
    /// Passing a single genome output prefix still works but is deprecated.
    #[arg(num_args = 0..=3, action = clap::ArgAction::Set, value_names = ["GENOME_SEQ", "GENOME_INS", "GENOME_DEL"])]
    pub genome_outputs: Vec<PathBuf>,

    /// The prefix to use for naming the output files (or an existing folder in
    /// which to place them).
    #[arg(long, conflicts_with = "sequence_output")]
    pub output_prefix: Option<PathBuf>,

    /// Name of the alignment module
    #[arg(short, long, default_value = "flu")]
    pub module: String,

    /// Run in simultaneous multi-threaded mode.
    #[arg(short = 'T', long)]
    pub threads: Option<NonZero<usize>>,

    /// Automatically detect the array task id from SGE or Slurm environment
    /// variables and write partition files for downstream collation.
    ///
    /// Output files are required and will be suffixed with a partition id.
    #[arg(short = 'G', long, conflicts_with_all = ["threads", "submit_grid_job"])]
    pub is_grid_task: bool,

    /// Submit and block on a grid engine (SGE or Slurm) array job of the
    /// specified size.
    #[arg(short = 'S', long, conflicts_with_all = ["threads", "is_grid_task"])]
    pub submit_grid_job: Option<usize>,

    /// Prints warning messages to stderr
    #[arg(long, conflicts_with_all = ["is_grid_task", "submit_grid_job"])]
    pub verbose: bool,

    /// A default ctype to use if any input records are not annotated. If not
    /// specified, an SSWSort module will be used to classify the query if a
    /// module exists, otherwise an error is produced.
    #[arg(long)]
    pub assume_default_ctype: Option<String>,
}

impl GridCompatibleCli for Cli {
    #[inline]
    fn is_grid_task(&self) -> bool {
        self.is_grid_task
    }

    #[inline]
    fn submit_grid_job(&self) -> Option<usize> {
        self.submit_grid_job
    }
}

/// A helper enum for differentiating the two major cases for the output files:
/// positional arguments are specified, or no positional arguments are specified
/// (in which case a provided or randomized prefix must be used).
enum CliOutputs {
    Positional {
        sequence_output:  PathBuf,
        insertion_output: PathBuf,
        deletion_output:  PathBuf,
        genome_outputs:   Vec<PathBuf>,
    },
    Prefix(PathBuf),
}

impl CliOutputs {
    /// Parses the outputs from `cli` into [`CliOutputs`], generating a random
    /// prefix and updating `matches` with it if needed.
    fn new(cli: &Cli, matches: &mut Vec<clap::ArgMatches>) -> Self {
        if let Some(sequence_output) = &cli.sequence_output {
            Self::Positional {
                sequence_output:  sequence_output.clone(),
                insertion_output: cli
                    .insertion_output
                    .as_ref()
                    .expect("Sequence output requires insertion output")
                    .clone(),
                deletion_output:  cli
                    .deletion_output
                    .as_ref()
                    .expect("Sequence output requires deletion output")
                    .clone(),
                genome_outputs:   cli.genome_outputs.clone(),
            }
        } else if let Some(prefix) = &cli.output_prefix {
            // Check whether the specified prefix is a directory, in which case
            // we must randomize the prefix within it
            if prefix.is_dir() {
                // Join output directory with random prefix
                let prefix = prefix.join(temp_name());

                // Clear previous value of --output-prefix
                for matches in matches.iter_mut() {
                    let _ = matches.try_clear_id("output_prefix");
                }

                // Inject the prefix as --output-prefix into the arg matches for
                // grid submission
                let temp_cmd = clap::Command::new("temp").arg(Arg::new("output_prefix").long("output-prefix").num_args(1));
                matches.push(temp_cmd.get_matches_from([
                    &OsStr::from("temp"),
                    &OsStr::from("--output-prefix"),
                    prefix.as_os_str(),
                ]));

                Self::Prefix(prefix)
            } else {
                Self::Prefix(prefix.clone())
            }
        } else {
            // Generate random prefix
            let prefix = temp_name();

            // Inject the prefix as --output-prefix into the arg matches for
            // grid submission
            let temp_cmd = clap::Command::new("temp").arg(Arg::new("output_prefix").long("output-prefix").num_args(1));
            matches.push(temp_cmd.get_matches_from(["temp", "--output-prefix", &prefix]));

            // Set output_prefix in Cli
            Self::Prefix(PathBuf::from(prefix))
        }
    }
}

/// A parsed version of [`Cli`] with all the output paths fully resolved.
pub struct Args {
    pub data_file:            PathBuf,
    pub module:               String,
    pub threads:              Option<NonZero<usize>>,
    pub submit_grid_job:      Option<usize>,
    pub verbose:              bool,
    pub assume_default_ctype: Option<String>,
    /// The sequence, insertion, and deletion output paths for the products
    pub product_output:       [PathBuf; 3],
    /// The sequence, insertion, and deletion output paths for the genome
    pub genome_output:        Option<[PathBuf; 3]>,
    /// The path to the log file to use for grid jobs
    pub grid_log_file:        PathBuf,
}

impl GridCompatibleArgs for Args {
    type Cli = Cli;

    fn from_cli(cli: Self::Cli) -> std::io::Result<Self> {
        Self::from_cli_with_matches(cli, &mut Vec::new())
    }

    fn from_cli_with_matches(cli: Cli, matches: &mut Vec<clap::ArgMatches>) -> std::io::Result<Self> {
        // Parse the outputs into one of two cases: positional arguments, or no
        // positional arguments and a prefix used instead
        let outputs = CliOutputs::new(&cli, matches);

        // Get the product outputs
        let product_output = match &outputs {
            CliOutputs::Positional {
                sequence_output,
                insertion_output,
                deletion_output,
                ..
            } => [sequence_output.clone(), insertion_output.clone(), deletion_output.clone()],
            CliOutputs::Prefix(path) => [
                path.with_added_extension("seq.txt"),
                path.with_added_extension("ins.txt"),
                path.with_added_extension("del.txt"),
            ],
        };

        // Get the genome outputs
        let genome_output = match &outputs {
            CliOutputs::Positional { genome_outputs, .. } => match genome_outputs.as_slice() {
                [] => None,
                [genome_prefix] => {
                    time_stamp(
                        "Warning: using the genome output prefix is deprecated, provide explicit genome output paths instead.",
                        true,
                    );
                    Some([
                        genome_prefix.to_owned(),
                        genome_prefix.with_added_extension("ins"),
                        genome_prefix.with_added_extension("del"),
                    ])
                }
                [seq, ins, del] => Some([seq.clone(), ins.clone(), del.clone()]),
                _ => {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "Expected either no genome output arguments, one genomic output prefix, or three genome output files.",
                    ));
                }
            },
            CliOutputs::Prefix(path) => Some([
                path.with_added_extension("gen_seq.txt"),
                path.with_added_extension("gen_ins.txt"),
                path.with_added_extension("gen_del.txt"),
            ]),
        };

        // For grid code, get the prefix used for the log file
        let grid_prefix = match outputs {
            CliOutputs::Positional { sequence_output, .. } => get_prefix(&sequence_output),
            CliOutputs::Prefix(path) => path,
        };

        // Get the log file path
        let grid_log_file = {
            let mut filename = grid_prefix
                .file_name()
                .ok_or(std::io::Error::other(format!(
                    "Failed to find filename in path: {}",
                    grid_prefix.display()
                )))?
                .to_os_string();
            filename.push("_ribosome_log.txt");
            grid_prefix.with_file_name(filename)
        };

        // Look for path duplicates in the upper triangle of paths
        let mut c1 = product_output.iter().chain(genome_output.iter().flatten());
        while let Some(p1) = c1.next() {
            let c2 = c1.clone();
            for p2 in c2 {
                if p1 == p2 {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "No two paths can be the same.",
                    ));
                }
            }
        }

        Ok(Self {
            data_file: cli.data_file,
            module: cli.module,
            threads: cli.threads,
            submit_grid_job: cli.submit_grid_job,
            verbose: cli.verbose,
            assume_default_ctype: cli.assume_default_ctype,
            product_output,
            genome_output,
            grid_log_file,
        })
    }

    fn outputs_mut(&mut self) -> impl Iterator<Item = &mut PathBuf> {
        self.product_output.iter_mut().chain(self.genome_output.iter_mut().flatten())
    }

    /// See [`GridCompatibleArgs::log_path`].
    ///
    /// ## Errors
    ///
    /// The `output_prefix` cannot end in `..`.
    fn log_path(&self) -> std::io::Result<PathBuf> {
        Ok(self.grid_log_file.clone())
    }
}

/// Helper to get a `mktemp`-like suffix for output
fn temp_name() -> String {
    let alpha = b"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
    let now = jiff::Timestamp::now().as_microsecond() as u64;
    let seed = getrandom::u64().unwrap_or(now);
    let seq = rand_sequence(alpha, 32, seed);
    String::from_utf8_lossy_owned(seq)
}

/// Infers the output prefix to use from the sequence file by stripping one
/// extension (followed by another if it is `seq`).
fn get_prefix(sequence_output: &Path) -> PathBuf {
    if sequence_output.extension().is_some() {
        let out = sequence_output.with_extension("");

        if out.extension().map(std::ffi::OsStr::as_encoded_bytes) == Some(b"seq") {
            out.with_extension("")
        } else {
            out
        }
    } else {
        sequence_output.to_path_buf()
    }
}
