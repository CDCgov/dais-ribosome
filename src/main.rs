#![feature(string_from_utf8_lossy_owned, bufreader_peek, try_trait_v2)]

use app::{
    args::Args,
    grid::{self, Grid},
    input::QueryInput,
    io, log,
};
use clap::Parser;
use dais_ribosome::{
    annotation::AnnotationModule,
    config::find_modules_toml,
    data::{ModuleData, RibosomeError},
};
use rayon::{iter::ParallelBridge, prelude::ParallelIterator};
use std::{collections::HashSet, path::Path};
use zoe::{data::err::OrFail, iter_utils::ProcessResultsExt};

// If we later want to match the shell script, we can use:
// <https://docs.rs/git-version/latest/git_version/macro.git_describe.html>
// const PROGRAM_VERSION: &str = concat!(env!("CARGO_PKG_NAME"), " v", env!("CARGO_PKG_VERSION"));

fn main() {
    let args = Args::parse();

    // Find the full file-system path to ribosome_res/modules.toml
    let toml_path = find_modules_toml().unwrap_or_fail();

    let mut module_data = ModuleData::load_from_file(&toml_path, &args.module)
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
        let g = Grid::task_vars_from_env();
        let writers = args.get_grid_writers(g.task_id).unwrap_or_fail();
        let gen_writers = args.get_optional_writers().unwrap_or_fail();

        process_queries_grid(&args.data_file, &annotation_module, writers, gen_writers, &g)
            .unwrap_or_die(&format!("Grid task {id} processing failed", id = g.task_id));
        return;
    }

    // Local execution
    let writers = args.get_writers().unwrap_or_fail();
    let gen_writers = args.get_optional_writers().unwrap_or_fail();

    let t = app::num_cpus::init_thread_pool(args.threads);
    if t > 1 {
        process_queries_parallel(&args.data_file, &annotation_module, writers, gen_writers)
            .unwrap_or_die("Query processing failed");
    } else {
        // Single-threaded, do not use a pool
        process_queries(&args.data_file, &annotation_module, writers, gen_writers).unwrap_or_die("Query processing failed");
    }
}

/// Processes queries sequentially (single-threaded).
fn process_queries(
    path: &Path, annotation_module: &AnnotationModule, mut writers: io::Writers, mut gen_writers: Option<io::Writers>,
) -> Result<(), RibosomeError> {
    log::ts("started, processing data");
    let queries = QueryInput::open(path)?;

    let mut unimplemented_ctypes = HashSet::new();

    for result in queries {
        let record = result?;
        let output = match annotation_module.process(record) {
            Ok(output) => output,
            Err(RibosomeError::UnimplementedCtype(ctype)) => {
                unimplemented_ctypes.insert(ctype);
                continue;
            }
            Err(e) => return Err(e),
        };

        let computed_output = output.materialize();

        io::write_output(&computed_output, &mut writers, &mut gen_writers)?;
    }

    log::print_unimplemented_ctypes(unimplemented_ctypes, annotation_module);
    log::ts("finished");

    Ok(())
}

/// Process queries in parallel using Rayon then write results sequentially.
fn process_queries_parallel(
    path: &Path, annotation_module: &AnnotationModule, mut writers: io::Writers, mut gen_writers: Option<io::Writers>,
) -> Result<(), RibosomeError> {
    log::ts("started, processing data (parallel)");
    let queries = QueryInput::open(path)?;

    let results: Vec<_> = queries.process_results(|iter| {
        iter.par_bridge()
            .map(|record| annotation_module.process(record).map(|o| o.materialize()))
            .collect::<Vec<_>>()
    })?;

    let mut unimplemented_ctypes = HashSet::new();

    for result in results {
        let computed_output = match result {
            Ok(output) => output,
            Err(RibosomeError::UnimplementedCtype(ctype)) => {
                unimplemented_ctypes.insert(ctype);
                continue;
            }
            Err(e) => return Err(e),
        };

        io::write_output(&computed_output, &mut writers, &mut gen_writers)?;
    }

    log::print_unimplemented_ctypes(unimplemented_ctypes, annotation_module);
    log::ts("finished");

    Ok(())
}

/// Process a grid partition: skip/step through queries based on task geometry.
fn process_queries_grid(
    path: &Path, annotation_module: &AnnotationModule, mut writers: io::Writers, mut gen_writers: Option<io::Writers>,
    g: &Grid,
) -> Result<(), RibosomeError> {
    let queries = QueryInput::open(path)?;

    // 0-based offset so we start on our partition
    let offset = (g.task_id - g.task_first) / g.task_stepsize;
    // modulus for interleaved partitioning
    let array_size = (g.task_last - g.task_first + 1).div_ceil(g.task_stepsize);

    for result in queries.skip(offset).step_by(array_size) {
        let record = result?;
        let output = match annotation_module.process(record) {
            Ok(output) => output,
            Err(RibosomeError::UnimplementedCtype(_)) => continue,
            Err(e) => return Err(e),
        };

        let computed_output = output.materialize();

        io::write_output(&computed_output, &mut writers, &mut gen_writers)?;
    }

    Ok(())
}

pub mod app;
