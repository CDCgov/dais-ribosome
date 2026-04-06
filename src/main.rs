#![feature(string_from_utf8_lossy_owned, bufreader_peek, try_trait_v2)]

use crate::app::{input::HandleNoNucleotidesExt, num_cpus::init_thread_pool, paths::find_modules_toml};
use app::{
    args::Args,
    grid::{self, Grid},
    input::QueryInput,
    log,
};
use clap::Parser;
use dais_ribosome::{
    AnnotationModule, ModuleData,
    error::{RibosomeError, UnimplementedCtype},
    toml::TomlConfig,
    tsv::{Writers, write_genome_output, write_product_output},
};
use rayon::{iter::ParallelBridge, prelude::ParallelIterator};
use std::{collections::HashSet, io::Write, path::Path};
use zoe::{data::err::OrFail, iter_utils::ProcessResultsExt};

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

    let annotation_module = module_data
        .build_annotation_module()
        .unwrap_or_die(&format!("Failed to build module '{}'", args.module));

    // Grid submission, blocks on the array
    if let Some(n) = args.submit_grid_job {
        let output_paths = args.output_paths_for_grid();
        grid::submit_job_sync(n, &args.module, &args.data_file, output_paths).unwrap_or_die("Grid job submission failed!");
        return;
    }

    // Grid array task
    if args.is_grid_task {
        let g = Grid::task_vars_from_env().unwrap_or_fail();
        let writers = args.get_grid_writers(g.task_id).unwrap_or_fail();
        let gen_writers = args.get_optional_writers().unwrap_or_fail();

        process_queries_grid(
            &args.data_file,
            &annotation_module,
            writers,
            gen_writers,
            &g,
            args.warn_no_nucleotides,
        )
        .unwrap_or_die(&format!("Grid task {id} processing failed", id = g.task_id));
        return;
    }

    // Local execution
    let writers = args.get_writers().unwrap_or_fail();
    let gen_writers = args.get_optional_writers().unwrap_or_fail();

    let num_threads = init_thread_pool(args.threads);

    if num_threads > 1 {
        process_queries_parallel(
            &args.data_file,
            &annotation_module,
            writers,
            gen_writers,
            args.warn_no_nucleotides,
        )
        .unwrap_or_die("Query processing failed");
    } else {
        // Single-threaded, do not use a pool
        process_queries(
            &args.data_file,
            &annotation_module,
            writers,
            gen_writers,
            args.warn_no_nucleotides,
        )
        .unwrap_or_die("Query processing failed");
    }
}

/// Processes queries sequentially (single-threaded).
fn process_queries<W: Write>(
    path: &Path, annotation_module: &AnnotationModule, mut writers: Writers<W>, mut gen_writers: Option<Writers<W>>,
    warn_no_nucleotides: bool,
) -> Result<(), RibosomeError> {
    log::ts("started, processing data");

    // Open the iterator of queries, and then remove any with an empty sequence
    // post-filtering, possibly issuing a warning
    let queries = QueryInput::open(path)?.handle_no_nucleotides(warn_no_nucleotides);

    let mut unimplemented_ctypes = HashSet::new();

    for result in queries {
        let record = result?;
        let output = match annotation_module.process(record) {
            Ok(output) => output,
            Err(RibosomeError::UnimplementedCtype(ctype)) => {
                unimplemented_ctypes.insert(ctype.0);
                continue;
            }
            Err(e) => return Err(e),
        };

        write_product_output(&output, &mut writers)?;
        if let Some(gen_writers) = &mut gen_writers {
            write_genome_output(&output, gen_writers)?;
        }
    }

    log::print_unimplemented_ctypes(unimplemented_ctypes, annotation_module);
    log::ts("finished");

    Ok(())
}

/// Process queries in parallel using Rayon then write results sequentially.
fn process_queries_parallel<W: Write>(
    path: &Path, annotation_module: &AnnotationModule, mut writers: Writers<W>, mut gen_writers: Option<Writers<W>>,
    warn_no_nucleotides: bool,
) -> Result<(), RibosomeError> {
    log::ts("started, processing data (parallel)");

    // Open the iterator of queries, and then remove any with an empty sequence
    // post-filtering, possibly issuing a warning. This is performed before
    // parallelization to avoid interleaved writes.
    let queries = QueryInput::open(path)?.handle_no_nucleotides(warn_no_nucleotides);

    // Use process_results to properly propagate catch errors that occur in the
    // input iterator. Collect into a result that propagates all errors instead
    // of UnimplementedCtype, so that rayon can end threads early when an error
    // occurs.

    let results = queries.process_results(|iter| {
        iter.par_bridge()
            .map(|record| match annotation_module.process(record) {
                Ok(output) => Ok(Ok(output)),
                Err(RibosomeError::UnimplementedCtype(e)) => Ok(Err(e)),
                Err(e) => Err(e),
            })
            .collect::<Result<Vec<_>, _>>()
    })??;

    let mut unimplemented_ctypes = HashSet::new();

    for result in results {
        let output = match result {
            Ok(output) => output,
            Err(UnimplementedCtype(ctype)) => {
                unimplemented_ctypes.insert(ctype);
                continue;
            }
        };

        write_product_output(&output, &mut writers)?;
        if let Some(gen_writers) = &mut gen_writers {
            write_genome_output(&output, gen_writers)?;
        }
    }

    log::print_unimplemented_ctypes(unimplemented_ctypes, annotation_module);
    log::ts("finished");

    Ok(())
}

/// Processes a grid partition: skip/step through queries based on task
/// geometry.
///
/// Note that this will not provide any messages regarding unimplemented
/// compound types.
fn process_queries_grid<W: Write>(
    path: &Path, annotation_module: &AnnotationModule, mut writers: Writers<W>, mut gen_writers: Option<Writers<W>>,
    g: &Grid, warn_no_nucleotides: bool,
) -> Result<(), RibosomeError> {
    // Open the iterator of queries, and then remove any with an empty sequence
    // post-filtering, possibly issuing a warning. This happens before
    // partitioning, so it may cause an offset between the n-th record in the
    // input file and the n-th record that gets partitioned.
    let queries = QueryInput::open(path)?.handle_no_nucleotides(warn_no_nucleotides);

    // 0-based offset so we start on our partition
    let offset = (g.task_id - g.task_first) / g.task_stepsize;
    // modulus for interleaved partitioning
    let array_size = (g.task_last - g.task_first + 1).div_ceil(g.task_stepsize);

    queries.process_results(|queries| {
        for record in queries.skip(offset).step_by(array_size) {
            let output = match annotation_module.process(record) {
                Ok(output) => output,
                Err(RibosomeError::UnimplementedCtype(_)) => continue,
                Err(e) => return Err(e),
            };

            write_product_output(&output, &mut writers)?;
            if let Some(gen_writers) = &mut gen_writers {
                write_genome_output(&output, gen_writers)?;
            }
        }

        Ok(())
    })?
}
