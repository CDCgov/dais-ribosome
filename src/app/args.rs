use crate::app::par_utils::grid::{GridCompatibleArgs, GridCompatibleCli};
use clap::{Arg, Parser};
use std::{num::NonZero, path::PathBuf};
use zoe::prelude::rand_sequence;

#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
/// CDS and and amino acid annotation tool for viruses.
pub struct Cli {
    /// Data file to annotate in TSV or FASTA format.
    ///         `†if classified, FASTA:  >ID|type_segment[_subtype]`
    ///         `*if classified, TSV:    ID<TAB>type_segment_[subtype]<TAB>sequence`
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

    /// Genomic file output prefix for sequences, insertion, and deletion.
    pub genomic_output_prefix: Option<PathBuf>,

    /// The prefix to use for naming the output files (or an existing folder in
    /// which to place them).
    #[arg(long, conflicts_with = "sequence_output")]
    pub output_prefix: Option<PathBuf>,

    /// Skips generated genome output when `--output-prefix` is specified.
    #[arg(long, conflicts_with = "genomic_output_prefix")]
    pub skip_genome: bool,

    /// Name of the alignment module
    #[arg(short, long, default_value = "flu")]
    pub module: String,

    /// Run in simultaneous multi-threaded mode.
    #[arg(short = 'T', long)]
    pub threads: Option<NonZero<usize>>,

    // TODO: This does nothing
    /// Write data as parquet files.
    #[arg(short = 'q', long)]
    pub output_parquet: bool,

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

    /// Prints warning messages to stderr. See the TODO for a full list of
    /// warnings that may be generated.
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

/// A parsed version of [`Cli`] with all the output paths fully resolved.
pub struct Args {
    pub data_file:            PathBuf,
    pub output_prefix:        PathBuf,
    pub module:               String,
    pub threads:              Option<NonZero<usize>>,
    pub output_parquet:       bool,
    pub submit_grid_job:      Option<usize>,
    pub verbose:              bool,
    pub assume_default_ctype: Option<String>,
    /// The sequence, insertion, and deletion output paths for the products
    pub product_output:       (PathBuf, PathBuf, PathBuf),
    /// The sequence, insertion, and deletion output paths for the genome
    pub genome_output:        Option<(PathBuf, PathBuf, PathBuf)>,
}

impl GridCompatibleArgs for Args {
    type Cli = Cli;

    fn from_cli(cli: Self::Cli) -> std::io::Result<Self> {
        Self::from_cli_with_matches(cli, &mut Vec::new())
    }

    fn from_cli_with_matches(cli: Cli, matches: &mut Vec<clap::ArgMatches>) -> std::io::Result<Self> {
        let output_prefix_or_dir = cli
            .output_prefix
            .into_iter()
            .chain(cli.genomic_output_prefix.clone())
            .chain(cli.sequence_output.as_ref().map(|p| p.with_extension("")))
            .chain(cli.insertion_output.as_ref().map(|p| p.with_extension("")))
            .chain(cli.deletion_output.as_ref().map(|p| p.with_extension("")))
            .next();

        let output_prefix = match output_prefix_or_dir {
            Some(output_dir) if output_dir.is_dir() => output_dir.join(temp_name()),
            Some(output_prefix) => output_prefix,
            None => {
                let prefix = temp_name();

                // Inject the prefix as --output-prefix into the arg matches for
                // grid submission
                let temp_cmd = clap::Command::new("temp").arg(Arg::new("output_prefix").long("output-prefix").num_args(1));
                matches.push(temp_cmd.get_matches_from(["temp", "--output-prefix", &prefix]));

                prefix.into()
            }
        };

        let sequence_output = cli
            .sequence_output
            .unwrap_or_else(|| output_prefix.with_added_extension("seq"));
        let insertion_output = cli
            .insertion_output
            .unwrap_or_else(|| output_prefix.with_added_extension("ins"));
        let deletion_output = cli
            .deletion_output
            .unwrap_or_else(|| output_prefix.with_added_extension("del"));

        let product_output = (sequence_output, insertion_output, deletion_output);

        let genome_prefix = cli
            .genomic_output_prefix
            .as_ref()
            .or((!cli.skip_genome).then_some(&output_prefix));

        let genome_output = if let Some(genome_prefix) = genome_prefix {
            let seq = genome_prefix.with_added_extension("gen_seq.txt");
            let ins = genome_prefix.with_added_extension("gen_ins.txt");
            let del = genome_prefix.with_added_extension("gen_del.txt");

            Some((seq, ins, del))
        } else {
            None
        };

        Ok(Self {
            data_file: cli.data_file,
            output_prefix,
            module: cli.module,
            threads: cli.threads,
            output_parquet: cli.output_parquet,
            submit_grid_job: cli.submit_grid_job,
            verbose: cli.verbose,
            assume_default_ctype: cli.assume_default_ctype,
            product_output,
            genome_output,
        })
    }

    fn outputs(&mut self) -> impl Iterator<Item = &mut PathBuf> {
        [
            &mut self.product_output.0,
            &mut self.product_output.1,
            &mut self.product_output.2,
        ]
        .into_iter()
        .chain(
            self.genome_output
                .iter_mut()
                .flat_map(move |genome_output| [&mut genome_output.0, &mut genome_output.1, &mut genome_output.2]),
        )
    }

    /// See [`GridCompatibleArgs::log_path`].
    ///
    /// ## Errors
    ///
    /// The `output_prefix` cannot end in `..`.
    fn log_path(&self) -> std::io::Result<PathBuf> {
        let mut file_name = self
            .output_prefix
            .file_name()
            .ok_or(std::io::Error::other(format!(
                "Failed to find filename in path: {}",
                self.output_prefix.display()
            )))?
            .to_os_string();
        file_name.push("_ribosome_log.txt");
        Ok(self.output_prefix.with_file_name(file_name))
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
