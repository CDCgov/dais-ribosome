//! Module data loading.

use crate::annotation::AnnotationModule;
use crate::config::toml::{AlignmentWeights, Formatting, Rules, TomlConfig};
use std::path::{Path, PathBuf};

use super::{
    ctype::build_ctype_map,
    error::ModuleLoadError,
    refs::{ReferenceMap, load_references},
    spec::{CdsSpecMap, load_cds_spec},
    weights::{CodonWeightMatrix, load_codon_weights},
};

// TODO: `cds_spec` is parsed into [`CdsSpecMap`], and then the first time it is
// used, it is converted into a different structure (map using RefKey instead).
// So, this could be altered.

/// Owned data backing an annotation module.
#[derive(Debug)]
pub struct ModuleData {
    pub name:              String,
    pub version:           String,
    pub formatting:        Formatting,
    pub rules:             Rules,
    pub alignment_weights: AlignmentWeights,
    pub other_modules:     Vec<(String, PathBuf)>,
    pub references:        ReferenceMap,
    // Consumed during build_annotation_module:
    cds_spec:              Option<CdsSpecMap>,
    codon_weights:         Option<CodonWeightMatrix>,
}

impl ModuleData {
    /// Load module data from a `modules.toml` configuration file.
    pub fn load_from_file(toml_path: &Path, module_name: &str) -> Result<Self, ModuleLoadError> {
        let config = TomlConfig::from_file(toml_path).map_err(|err| ModuleLoadError::io(toml_path, err))?;
        Self::from_config(toml_path, config, module_name)
    }

    /// Build module data from a parsed configuration.
    pub fn from_config(toml_path: &Path, config: TomlConfig, module_name: &str) -> Result<Self, ModuleLoadError> {
        // Get path to ribosome_res directory
        let modules_dir = toml_path
            .parent()
            .expect("The modules.toml path must have parent ribosome_res");

        // Select the module of interest given module_name
        let (module, other_modules) = config
            .find_module(module_name, modules_dir)
            .ok_or_else(|| ModuleLoadError::module_not_found(module_name))?;

        // Get path to module folder
        let module_root = modules_dir.join(&module.name);

        // Get paths to files within module folder
        let references_path = module_root.join(&module.references);
        let weights_path = module_root.join(&module.weights);
        let cds_spec_path = module_root.join(&module.cds_spec);

        let references = load_references(&references_path).map_err(|err| ModuleLoadError::io(&references_path, err))?;
        let codon_weights = load_codon_weights(&weights_path).map_err(|err| ModuleLoadError::io(&weights_path, err))?;
        let cds_spec = load_cds_spec(&cds_spec_path).map_err(|err| ModuleLoadError::io(&cds_spec_path, err))?;
        let alignment_weights = module.alignment;

        Ok(Self {
            name: module.name,
            version: module.version.unwrap_or_else(|| "unknown".to_string()),
            formatting: module.formatting,
            rules: module.rules,
            alignment_weights,
            other_modules,
            references,
            cds_spec: Some(cds_spec),
            codon_weights: Some(codon_weights),
        })
    }

    // TODO: Weird design, would be good to change if possible.

    /// Builds an [`AnnotationModule`] that borrows from this [`ModuleData`].
    /// This consumes `cds_spec` and `codon_weights`, so it can only be called
    /// once.
    pub fn build_annotation_module(&mut self) -> Result<AnnotationModule<'_>, ModuleLoadError> {
        let cds_spec = self
            .cds_spec
            .take()
            .ok_or_else(|| ModuleLoadError::validation("build_annotation_module can only be called once"))?;

        let codon_weights = self
            .codon_weights
            .take()
            .ok_or_else(|| ModuleLoadError::validation("build_annotation_module can only be called once"))?;
        let ctype_map = build_ctype_map(&self.references, cds_spec, codon_weights, &self.alignment_weights)?;
        Ok(AnnotationModule { data: self, ctype_map })
    }
}
