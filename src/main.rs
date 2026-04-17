#![feature(string_from_utf8_lossy_owned, bufreader_peek, try_trait_v2, iter_intersperse)]

use crate::app::{
    input::{NoCtype, QueryInfo},
    log::time_stamp,
    num_cpus::init_thread_pool,
    paths::find_modules_toml,
};
use app::{
    args::Args,
    grid::{self, Grid},
    input::QueryReader,
    log,
};
use clap::Parser;
use dais_ribosome::{
    AnnotationModule, ModuleData,
    error::{RibosomeError, UnimplementedCtype},
    outputs::RibosomeOutput,
    toml::TomlConfig,
    tsv::{Writers, write_genome_output, write_product_output},
};
use rayon::{iter::ParallelBridge, prelude::ParallelIterator};
use sswsort::SSWSortModule;
use std::{collections::HashSet, error::Error, fmt::Display, io::Write, path::PathBuf};
use zoe::{
    data::err::{Fail, GetCode, OrFail, ResultWithErrorContext, WithErrorContext},
    iter_utils::ProcessResultsExt,
};

pub mod app;

// If we later want to match the shell script, we can use:
// <https://docs.rs/git-version/latest/git_version/macro.git_describe.html>
// const PROGRAM_VERSION: &str = concat!(env!("CARGO_PKG_NAME"), " v", env!("CARGO_PKG_VERSION"));

fn main() {
    let args = Args::parse();

    // Find the full file-system path to ribosome_res/modules.toml
    let toml_path = find_modules_toml().unwrap_or_fail();

    // Parse the TOML file
    let parsed_toml = TomlConfig::from_file(&toml_path).unwrap_or_fail();

    // Convert the TOML data into ModuleData
    let module_data = ModuleData::new(parsed_toml, &toml_path, &args.module)
        .unwrap_or_die(&format!("Failed to prepare module '{}'", args.module));

    // Build the AnnotationModule
    let annotation_module = module_data
        .build_annotation_module()
        .unwrap_or_die(&format!("Failed to build module '{}'", args.module));

    let classification = ClassificationStrategy::new(&args).unwrap_or_fail();

    // Grid submission, blocks on the array
    if let Some(n) = args.submit_grid_job {
        let output_paths = args.output_paths_for_grid();
        grid::submit_job_sync(n, &args.module, &args.data_file, output_paths).unwrap_or_die("Grid job submission failed!");
        return;
    }

    let config = BinaryConfig {
        classification,
        annotation: annotation_module,
        verbose: args.verbose,
        list_unimplemented_ctypes: !args.is_grid_task,
    };

    let queries = QueryReader::from_path(&args.data_file).unwrap_or_fail();

    let result = if args.is_grid_task {
        let g = Grid::task_vars_from_env().unwrap_or_fail();
        let writers = args.get_grid_writers(g.task_id).unwrap_or_fail();
        let gen_writers = args.get_optional_writers().unwrap_or_fail();

        // 0-based offset so we start on our partition
        let offset = (g.task_id - g.task_first) / g.task_stepsize;
        // modulus for interleaved partitioning
        let array_size = (g.task_last - g.task_first + 1).div_ceil(g.task_stepsize);

        queries.process_results(|queries| {
            let stepped = queries.skip(offset).step_by(array_size);
            process_queries(stepped, &config, writers, gen_writers)
        })
    } else {
        let writers = args.get_writers().unwrap_or_fail();
        let gen_writers = args.get_optional_writers().unwrap_or_fail();

        let num_threads = init_thread_pool(args.threads);

        if num_threads > 1 {
            queries.process_results(|queries| process_queries_parallel(queries, &config, writers, gen_writers))
        } else {
            queries.process_results(|queries| process_queries(queries, &config, writers, gen_writers))
        }
    };

    match result {
        Ok(Ok(())) => {}
        Ok(Err(ProcessingError::NoCtype(e))) => {
            let err = std::io::Error::other(e.with_context(format!(
                "All ctypes must be specified since SSWSort does not have a module for: {module}",
                module = args.module
            )));
            err.fail()
        }
        Ok(Err(e)) => e.fail(),
        Err(e) => e.fail(),
    }
}

/// Configuration for the binary portion of DAIS-ribosome.
///
/// This groups together the modules for SSWSort, DAIS-ribosome, as well as
/// other configuration.
pub struct BinaryConfig<'a> {
    /// The strategy to use for assigning a `ctype` to unannotated inputs.
    classification:            Option<ClassificationStrategy>,
    /// The module to use for performing annotation and translation.
    annotation:                AnnotationModule<'a>,
    /// Whether to display warnings.
    verbose:                   bool,
    /// Whether to list any unimplemented compound types encountered (i.e.,
    /// those not handled by the annotation module).
    list_unimplemented_ctypes: bool,
}

/// The possible ways of handling a missing ctype in an incoming record.
pub enum ClassificationStrategy {
    /// Attempt to classify the ctype using the given SSWSort module.
    SswSort(SSWSortModule),
    /// Use a default ctype.
    Default(String),
}

impl ClassificationStrategy {
    /// Initializes the [`ClassificationStrategy`] based on the passed [`Args`].
    ///
    /// If `--assume-default-ctype` is passed, then
    /// [`ClassificationStrategy::Default`] is used. Otherwise, the function
    /// attempts to locate an SSWSort module of the same name as the
    /// DAIS-ribosome module. If that fails, `None` is returned.
    ///
    /// ## Errors
    ///
    /// If the SSWSort TOML path exists, then parsing or IO errors are
    /// propagated. If a module with the requested name is found, then any
    /// errors opening the references are also propagated. Context is added to
    /// all errors.
    pub fn new(args: &Args) -> std::io::Result<Option<Self>> {
        if let Some(default) = &args.assume_default_ctype {
            return Ok(Some(ClassificationStrategy::Default(default.clone())));
        }

        let sswsort_toml_path = PathBuf::from("sswsort_res/config.toml");

        let sswsort_module = if sswsort_toml_path.exists() {
            let config = sswsort::TomlConfig::from_path(&sswsort_toml_path)
                .with_path_context("Failed to parse SSWSort TOML", sswsort_toml_path)?;
            if let Some(params) = config.get(&args.module) {
                Some(SSWSortModule::new(params).with_context("Failed to load SSWSort reference sequences")?)
            } else {
                None
            }
        } else {
            None
        };

        Ok(sswsort_module.map(ClassificationStrategy::SswSort))
    }
}

/// An error that could arise from [`process_queries`] or
/// [`process_queries_parallel`].
///
/// This error type is specific to the binary.
#[derive(Debug)]
pub enum ProcessingError {
    /// A ctype was not specified in an input file, and no classification
    /// strategy was specified.
    NoCtype(NoCtype),
    /// A file was empty.
    EmptyFile(std::path::PathBuf),
    /// Any other error that may have occurred.
    Io(std::io::Error),
}

impl From<std::io::Error> for ProcessingError {
    fn from(value: std::io::Error) -> Self {
        ProcessingError::Io(value)
    }
}

impl From<NoCtype> for ProcessingError {
    fn from(value: NoCtype) -> Self {
        ProcessingError::NoCtype(value)
    }
}

impl Display for ProcessingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ProcessingError::NoCtype(e) => e.fmt(f),
            ProcessingError::EmptyFile(p) => write!(f, "Empty file: {}", p.display()),
            ProcessingError::Io(e) => e.fmt(f),
        }
    }
}

impl Error for ProcessingError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            ProcessingError::NoCtype(e) => e.source(),
            ProcessingError::EmptyFile(_) => None,
            ProcessingError::Io(e) => e.source(),
        }
    }
}

impl GetCode for ProcessingError {
    fn get_code(&self) -> i32 {
        match self {
            ProcessingError::NoCtype(e) => e.get_code(),
            ProcessingError::EmptyFile(_) => 66, // EX_NOINPUT
            ProcessingError::Io(e) => e.get_code(),
        }
    }
}

/// Processes queries sequentially (single-threaded).
fn process_queries<Q, W>(
    queries: Q, config: &BinaryConfig, mut writers: Writers<W>, mut gen_writers: Option<Writers<W>>,
) -> Result<(), ProcessingError>
where
    Q: Iterator<Item = QueryInfo>,
    W: Write, {
    log::ts("started, processing data");

    let mut unimplemented_ctypes = HashSet::new();

    for info in queries {
        let record = match info.classify_and_prepare(&config.classification, config.verbose) {
            Ok(Some(record)) => record,
            Ok(None) => continue,
            Err(e) => return Err(ProcessingError::NoCtype(e)),
        };

        let output = match config.annotation.process(record) {
            Ok(output) => output,
            Err(RibosomeError::UnimplementedCtype(ctype)) => {
                unimplemented_ctypes.insert(ctype.0);
                continue;
            }
            Err(RibosomeError::EmptyFile(e)) => return Err(ProcessingError::EmptyFile(e)),
            Err(RibosomeError::Io(e)) => return Err(ProcessingError::Io(e)),
        };

        if config.verbose {
            warn_failed_ref_ids(&output);
        }

        write_product_output(&output, &mut writers)?;
        if let Some(gen_writers) = &mut gen_writers {
            write_genome_output(&output, gen_writers)?;
        }
    }

    if config.list_unimplemented_ctypes {
        log::print_unimplemented_ctypes(unimplemented_ctypes, &config.annotation);
    }

    log::ts("finished");

    Ok(())
}

/// Process queries in parallel using Rayon then write results sequentially.
fn process_queries_parallel<Q, W>(
    queries: Q, config: &BinaryConfig, mut writers: Writers<W>, mut gen_writers: Option<Writers<W>>,
) -> Result<(), ProcessingError>
where
    Q: Iterator<Item = QueryInfo> + ParallelBridge + Send,
    W: Write, {
    log::ts("started, processing data (parallel)");

    // Use process_results to properly propagate catch errors that occur in the
    // input iterator. Collect into a result that propagates all errors instead
    // of UnimplementedCtype, so that rayon can end threads early when an error
    // occurs.

    let results = queries
        .par_bridge()
        .flat_map(|info| info.classify_and_prepare(&config.classification, config.verbose).transpose())
        .map(|record| match config.annotation.process(record?) {
            Ok(output) => Ok(Ok(output)),
            Err(RibosomeError::UnimplementedCtype(e)) => Ok(Err(e)),
            Err(RibosomeError::EmptyFile(e)) => Err(ProcessingError::EmptyFile(e)),
            Err(RibosomeError::Io(e)) => Err(ProcessingError::Io(e)),
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut unimplemented_ctypes = HashSet::new();

    for result in results {
        let output = match result {
            Ok(output) => output,
            Err(UnimplementedCtype(ctype)) => {
                unimplemented_ctypes.insert(ctype);
                continue;
            }
        };

        if config.verbose {
            warn_failed_ref_ids(&output);
        }

        write_product_output(&output, &mut writers)?;
        if let Some(gen_writers) = &mut gen_writers {
            write_genome_output(&output, gen_writers)?;
        }
    }

    if config.list_unimplemented_ctypes {
        log::print_unimplemented_ctypes(unimplemented_ctypes, &config.annotation);
    }
    log::ts("finished");

    Ok(())
}

fn warn_failed_ref_ids(output: &RibosomeOutput) {
    for ref_id in &output.failed_ref_ids {
        time_stamp(
            &format!(
                "Failed to align query {query_id} against reference {ref_id}",
                query_id = output.query.id(),
            ),
            true,
        );
    }
}
