//! Protein annotation engine.
pub mod error;
pub mod hashing;
mod process;

use crate::data::{ctype::CompoundTypeMap, keys::RefKey, module::ModuleData};
use std::collections::{HashMap, HashSet};
use zoe::prelude::*;

/// A fully materialized annotation module ready for protein annotation work.
#[derive(Debug)]
pub struct AnnotationModule<'a> {
    /// Reference to the backing module data (for weight matrices, etc.).
    pub(crate) data:      &'a ModuleData,
    /// Compound type map for iteration-based processing.
    pub(crate) ctype_map: CompoundTypeMap<'a>,
}

impl<'a> AnnotationModule<'a> {
    /// Suggest modules that might contain given compound types.
    pub fn suggest_modules_for_compound_types(&self, mut compound_types: HashSet<String>) -> HashMap<&str, Vec<String>> {
        let mut found: HashMap<&str, Vec<String>> = HashMap::new();
        for (module_name, ref_path) in &self.data.other_modules {
            if let Ok(reader) = FastaReader::from_filename(ref_path) {
                for fasta_result in reader.flatten() {
                    if let Some(key) = RefKey::parse(&fasta_result.name)
                        && let Some(ctype) = compound_types.take(&key.compound_type)
                    {
                        found.entry(module_name.as_str()).or_default().push(ctype);
                    }
                }
            }
        }
        found
    }
}
