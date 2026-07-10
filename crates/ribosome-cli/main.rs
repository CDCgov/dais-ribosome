#![feature(string_from_utf8_lossy_owned, bufreader_peek, try_trait_v2, iter_intersperse)]

use args::Args;
use dais_ribosome::{AnnotationModule, errors::RibosomeError, outputs::RibosomeOutput, toml::TomlConfig, tsv::Writers};
use input::{NoCtype, QueryInfo, QueryReader};
use log::time_stamp;
use num_cpus::init_thread_pool;
use par_utils::{
    grid::{GridCompatibleArgs, GridInfo, JobErrorOrFail},
    writers::WriterThreaded,
};
use paths::find_modules_toml;
use rayon::{iter::ParallelBridge, prelude::ParallelIterator};
use sswsort::SSWSortModule;
use std::{collections::HashSet, error::Error, fmt::Display, io::Write, path::PathBuf};
use zoe::{
    data::err::{Fail, GetCode, OrFail, ResultWithErrorContext, WithErrorContext},
    iter_utils::ProcessResultsExt,
    unwrap_or_return_some_err,
};

pub mod args;
pub mod input;
pub mod log;
pub mod num_cpus;
pub mod par_utils;
pub mod paths;

// If we later want to match the shell script, we can use:
// <https://docs.rs/git-version/latest/git_version/macro.git_describe.html>
// const PROGRAM_VERSION: &str = concat!(env!("CARGO_PKG_NAME"), " v", env!("CARGO_PKG_VERSION"));

fn main() {
    // Parse the arguments, get grid info, adjust paths based on task ID, open
    // writers.
    let (args, grid_info) = Args::parse_maybe_grid().unwrap_or_fail();

    // Find the full file-system path to ribosome_res/modules.toml
    let toml_path = find_modules_toml().unwrap_or_fail();

    // Parse the TOML file
    let parsed_toml = TomlConfig::from_file(&toml_path).unwrap_or_fail();

    // Build the AnnotationModule
    let annotation_module = AnnotationModule::new(&parsed_toml, &toml_path, &args.module)
        .unwrap_or_die(&format!("Failed to build module '{}'", args.module));

    // Determine the classification strategy, which may involve loading SSWSort
    // module
    let classification = ClassificationStrategy::new(&args).unwrap_or_fail();

    // Handle a request to submit a grid job
    let grid_info = match grid_info {
        Some(GridInfo::Requested(grid_info)) => {
            grid_info.submit_job_sync().unwrap_or_fail();
            return;
        }
        Some(GridInfo::Task(grid_info)) => Some(grid_info),
        None => None,
    };

    // Group relevant information together for ease of function calls
    let config = BinaryConfig {
        classification,
        annotation: annotation_module,
        verbose: args.verbose,
        list_unimplemented_ctypes: grid_info.is_none(),
    };

    // Initialize an iterator over the input queries
    let queries = QueryReader::from_path(&args.data_file)
        .with_path_context("Failed to open query file", args.data_file)
        .map_err(std::io::Error::from)
        .unwrap_or_fail();

    // Open the writers
    let [seq, ins, del] = args.product_output;
    let writers = Writers::from_paths(seq, ins, del).unwrap_or_die("Failed to create product output files");
    let gen_writers = if let Some([seq, ins, del]) = args.genome_output {
        Some(Writers::from_paths(seq, ins, del).unwrap_or_die("Failed to create genome output files"))
    } else {
        None
    };

    // Handle processing in the case of a grid task
    if let Some(grid_info) = grid_info {
        grid_info.run_task(|grid_info| {
            log::ts("started, processing data (grid)");

            let result = queries.process_results(|queries| {
                let stepped = grid_info.select_inputs(queries);
                process_queries(stepped, &config, writers, gen_writers)
            });

            // Handle any errors in the query iterator or the processing
            match result {
                Ok(Ok(())) => {}
                Ok(Err(ProcessingError::NoCtype(err))) => {
                    log::ts("annotation failed");
                    let err = std::io::Error::from(err.with_context(format!(
                        "All ctypes must be specified since SSWSort does not have a module for: {module}",
                        module = args.module
                    )));
                    err.fail()
                }
                Ok(Err(err)) => {
                    log::ts("annotation failed");
                    err.fail()
                }
                Err(err) => {
                    log::ts("annotation failed");
                    err.fail()
                }
            }

            log::ts("finished");
        });
    }

    // Get the number of threads
    let num_threads = init_thread_pool(args.threads);

    // Handle processing for either single threaded or multithreaded
    let result = if num_threads > 1 {
        log::ts("started, processing data (parallel)");
        queries.process_results(|queries| process_queries_parallel(queries, &config, writers, gen_writers))
    } else {
        log::ts("started, processing data");
        queries.process_results(|queries| process_queries(queries, &config, writers, gen_writers))
    };

    // Handle any errors in the query iterator or the processing
    match result {
        Ok(Ok(())) => {}
        Ok(Err(ProcessingError::NoCtype(err))) => {
            log::ts("annotation failed");
            let err = std::io::Error::from(err.with_context(format!(
                "All ctypes must be specified since SSWSort does not have a module for: {module}",
                module = args.module
            )));
            err.fail()
        }
        Ok(Err(err)) => {
            log::ts("annotation failed");
            err.fail()
        }
        Err(err) => {
            log::ts("annotation failed");
            err.fail()
        }
    }

    log::ts("finished");
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
            ProcessingError::Io(e) => e.fmt(f),
        }
    }
}

impl Error for ProcessingError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            ProcessingError::NoCtype(e) => e.source(),
            ProcessingError::Io(e) => e.source(),
        }
    }
}

impl GetCode for ProcessingError {
    fn get_code(&self) -> i32 {
        match self {
            ProcessingError::NoCtype(e) => e.get_code(),
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
                unimplemented_ctypes.insert(ctype);
                continue;
            }
            Err(RibosomeError::Io(e)) => return Err(ProcessingError::Io(e)),
        };

        if config.verbose {
            warn_failed_ref_ids(&output);
        }

        writers.write_product_output(&output)?;
        if let Some(gen_writers) = &mut gen_writers {
            gen_writers.write_genome_output(&output)?;
        }
    }

    writers.flush()?;
    if let Some(mut gen_writers) = gen_writers {
        gen_writers.flush()?;
    }

    if config.list_unimplemented_ctypes {
        log::print_unimplemented_ctypes(unimplemented_ctypes, &config.annotation);
    }

    Ok(())
}

/// Process queries in parallel using Rayon then write results sequentially.
fn process_queries_parallel<Q, W>(
    queries: Q, config: &BinaryConfig, writers: Writers<W>, gen_writers: Option<Writers<W>>,
) -> Result<(), ProcessingError>
where
    Q: Iterator<Item = QueryInfo> + ParallelBridge + Send,
    W: Write + Send + 'static, {
    // Use process_results to properly propagate catch errors that occur in the
    // input iterator. Collect into a result that propagates all errors instead
    // of UnimplementedCtype, so that rayon can end threads early when an error
    // occurs.

    let mut writers = writers.map(WriterThreaded::new);
    let gen_writers = gen_writers.map(|gen_writers| gen_writers.map(WriterThreaded::new));

    let unimplemented_ctypes = queries
        .par_bridge()
        .flat_map(|info| info.classify_and_prepare(&config.classification, config.verbose).transpose())
        .map_with((writers.clone(), gen_writers.clone()), |(writers, gen_writers), record| {
            let record = unwrap_or_return_some_err!(record.map_err(Into::into));

            match config.annotation.process(record) {
                Ok(output) => {
                    if config.verbose {
                        warn_failed_ref_ids(&output);
                    }

                    unwrap_or_return_some_err!(writers.write_product_output(&output).map_err(Into::into));
                    if let Some(gen_writers) = gen_writers {
                        unwrap_or_return_some_err!(gen_writers.write_genome_output(&output).map_err(Into::into));
                    }
                    None
                }
                Err(RibosomeError::UnimplementedCtype(e)) => Some(Ok(e)),
                Err(RibosomeError::Io(e)) => Some(Err(ProcessingError::Io(e))),
            }
        })
        .flatten()
        .collect::<Result<HashSet<_>, _>>()?;

    writers.flush()?;
    if let Some(mut gen_writers) = gen_writers {
        gen_writers.flush()?;
    }

    if config.list_unimplemented_ctypes {
        log::print_unimplemented_ctypes(unimplemented_ctypes, &config.annotation);
    }

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
