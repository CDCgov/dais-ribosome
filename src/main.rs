#![feature(string_from_utf8_lossy_owned, bufreader_peek, try_trait_v2)]
//use rayon::{ThreadPoolBuilder, iter::ParallelBridge, prelude::ParallelIterator};

use app::{args::Args, input::QueryInput, io::Writers, log};

use clap::Parser;
use dais_ribosome::{
    annotation::AnnotationModule,
    config::find_modules_toml,
    data::{ModuleData, RibosomeError},
};
use std::{collections::HashSet, io::Write, path::Path};
use zoe::data::err::OrFail;

// If we later want to match the shell script, we can use:
// <https://docs.rs/git-version/latest/git_version/macro.git_describe.html>
// const PROGRAM_VERSION: &str = concat!(env!("CARGO_PKG_NAME"), " v", env!("CARGO_PKG_VERSION"));

fn main() {
    let args = Args::parse();

    let toml_path = find_modules_toml().unwrap_or_fail();

    let mut module_data = ModuleData::load_from_file(&toml_path, &args.module)
        .unwrap_or_die(&format!("Failed to prepare module '{}'", args.module));

    let annotation_module = module_data
        .build_annotation_module()
        .unwrap_or_die(&format!("Failed to build module '{}'", args.module));

    let writers = args.get_writers().unwrap_or_fail();
    let gen_writers = args.get_optional_writers().unwrap_or_fail();

    process_queries(&args.data_file, &annotation_module, writers, gen_writers).unwrap_or_die("Query processing failed");
}

// TODO: move later
fn process_queries(
    path: &Path, annotation_module: &AnnotationModule, mut writers: Writers, mut gen_writers: Option<Writers>,
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

        for row in output.seq_rows() {
            writeln!(writers.seq, "{row}")?;
        }

        for row in output.ins_rows() {
            writeln!(writers.ins, "{row}")?;
        }

        for row in output.del_rows() {
            writeln!(writers.del, "{row}")?;
        }

        // Genome output rows
        if let Some(ref mut w) = gen_writers {
            for row in output.gen_rows() {
                writeln!(w.seq, "{row}")?;
            }

            for row in output.gen_ins_rows() {
                writeln!(w.ins, "{row}")?;
            }

            for row in output.gen_del_rows() {
                writeln!(w.del, "{row}")?;
            }
        }
    }

    log::print_unimplemented_ctypes(unimplemented_ctypes, annotation_module);
    log::ts("finished");

    Ok(())
}

pub mod app;
