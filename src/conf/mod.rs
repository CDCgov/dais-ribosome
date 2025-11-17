mod refs;
mod spec;
mod toml;
mod weights;

pub use refs::*;
pub use spec::*;
pub use toml::*;
pub use weights::*;

use zoe::data::records::RecordReader;

pub struct TSVReader {}
impl RecordReader for TSVReader {
    const RECORD_NAME: &str = "TSV";
}

pub struct CompoundType(pub String);
