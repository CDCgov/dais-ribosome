//! Module data loading.

use crate::{
    annotation::AnnotationModule,
    config::{AlignmentWeights, Formatting, Rules, TomlConfig},
    data::{
        ctype::build_ctype_map,
        refs::{ReferenceMap, load_references},
        spec::{CdsSpecMap, load_cds_spec},
        weights::{CodonWeightMatrix, load_codon_weights},
    },
};
use std::path::{Path, PathBuf};
use zoe::data::err::ResultWithErrorContext;

/// Owned data backing an annotation module.
#[derive(Debug)]
pub struct ModuleData {
    /// The name of the module (e.g., `flu`, `cov`, or `rsv`). This must
    /// correspond to a folder in `ribosome_res`.
    pub name:              String,
    /// An optional version for the module (e.g., `2.0-alpha`).
    pub version:           String,
    pub formatting:        Formatting,
    pub rules:             Rules,
    /// Alignments weights for the module.
    pub alignment_weights: AlignmentWeights,
    /// A vector of the other module names and the paths to their reference
    /// files, used in for providing warning messages when unrecognized compound
    /// types are encountered.
    pub other_modules:     Vec<(String, PathBuf)>,
    /// A hash map from [`RefKey`] values to a vector of reference sequences.
    ///
    /// [`RefKey`]: crate::data::keys::RefKey
    pub references:        ReferenceMap,

    // Consumed during build_annotation_module:
    /// A hash map from [`RefKey`] values to a vector of the protein product
    /// names (e.g., `HA`, `HA-signal`) and their [`Exons`].
    ///
    /// [`RefKey`]: crate::data::keys::RefKey
    /// [`Exons`]: crate::data::exons::Exons
    cds_spec:      Option<CdsSpecMap>,
    codon_weights: Option<CodonWeightMatrix>,
}

impl ModuleData {
    /// Build module data from a parsed configuration.
    ///
    /// ## Errors
    ///
    /// Any errors are returned without including the `module_name` as context,
    /// since the caller will add that. The `toml_path` is included as context
    /// when relevant.
    pub fn new(config: TomlConfig, toml_path: &Path, module_name: &str) -> std::io::Result<Self> {
        // Get path to ribosome_res directory
        let modules_dir = toml_path
            .parent()
            .expect("The modules.toml path must have parent ribosome_res");

        // Select the module of interest given module_name
        let Some((module, other_modules)) = config.find_module(module_name, modules_dir) else {
            return Err(std::io::Error::other(format!(
                "Module name not found in configuration file: {toml_path}",
                toml_path = toml_path.display()
            )));
        };

        // Get path to module folder
        let module_root = modules_dir.join(&module.name);

        // Get paths to files within module folder
        let references_path = module_root.join(&module.references);
        let weights_path = module_root.join(&module.weights);
        let cds_spec_path = module_root.join(&module.cds_spec);

        let references = load_references(&references_path)
            .with_path_context("Failed to load the references from file", &references_path)?;
        let codon_weights = load_codon_weights(&weights_path)
            .with_path_context("Failed to load the codon position weights from file", weights_path)?;
        let cds_spec =
            load_cds_spec(&cds_spec_path).with_path_context("Failed to load the CDS specs from file", cds_spec_path)?;
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
    pub fn build_annotation_module(&mut self) -> std::io::Result<AnnotationModule<'_>> {
        let cds_spec = self.cds_spec.take().expect("build_annotation_module can only be called once");

        let codon_weights = self
            .codon_weights
            .take()
            .expect("build_annotation_module can only be called once");

        let ctype_map = build_ctype_map(&self.references, cds_spec, codon_weights, &self.alignment_weights)?;
        Ok(AnnotationModule { data: self, ctype_map })
    }
}
