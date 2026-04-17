use crate::app::grid::get_partition_filename;
use dais_ribosome::{error::RibosomeError, tsv::Writers};
use std::{fs::File, io::BufWriter, path::PathBuf, sync::LazyLock};
use zoe::{data::err::ResultWithErrorContext, prelude::rand_sequence};

pub fn open_writer(path: &Option<PathBuf>, extension: &str) -> Result<BufWriter<File>, RibosomeError> {
    let mut new_path = PathBuf::with_capacity(TEMP_PREFIX.len() + extension.len() + 1);
    let p = path.as_ref().unwrap_or_else(|| {
        new_path.push(&*TEMP_PREFIX);
        new_path.add_extension(extension);
        &new_path
    });
    let f = File::create(p).with_path_context(format!("Could not open writer for {extension} file"), p)?;
    Ok(BufWriter::new(f))
}

/// Open a partition writer for grid array tasks.
///
/// The partition filename is derived from the output path with a `_NNN` suffix.
pub fn open_partition_writer(
    path: &Option<PathBuf>, task_id: usize, extension: &str,
) -> Result<BufWriter<File>, RibosomeError> {
    let p = path.as_ref().expect("Grid task mode requires explicit output paths");
    let partition_name = get_partition_filename(p, task_id, extension);
    let mut partition_path = p.clone();
    partition_path.set_file_name(partition_name);
    let f = File::create(&partition_path).with_path_context(
        format!("Could not open partition writer for {extension} file"),
        &partition_path,
    )?;
    Ok(BufWriter::new(f))
}

pub fn optional_writers(path: &Option<PathBuf>) -> Result<Option<Writers<BufWriter<File>>>, RibosomeError> {
    if let Some(p) = path {
        let mut p = p.clone();
        p.set_extension("gen_seq.txt");
        let seq = BufWriter::new(File::create(&p).with_path_context("Could not open writer for genome sequence file", &p)?);

        p.set_extension("");
        p.set_extension("");
        p.set_extension("gen_ins.txt");
        let ins = BufWriter::new(File::create(&p).with_path_context(
            format!("Could not open writer for genome insertion file: {}", p.display()),
            &p,
        )?);

        p.set_extension("");
        p.set_extension("");
        p.set_extension("gen_del.txt");
        let del = BufWriter::new(File::create(&p).with_path_context("Could not open writer for genome deletion file", &p)?);

        Ok(Some(Writers { seq, ins, del }))
    } else {
        Ok(None)
    }
}

static TEMP_PREFIX: LazyLock<String> = LazyLock::new(temp_name);

/// Helper to get a `mktemp`-like suffix for output
fn temp_name() -> String {
    let alpha = b"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
    let now = jiff::Timestamp::now().as_microsecond() as u64;
    let seed = getrandom::u64().unwrap_or(now);
    let seq = rand_sequence(alpha, 32, seed);
    String::from_utf8_lossy_owned(seq)
}
