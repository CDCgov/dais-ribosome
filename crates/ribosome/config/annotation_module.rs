//! The [`AnnotationModule`] for containing all data needed to translate and
//! annotate query sequences, along with helper structs.

use crate::{
    config::{
        ProductSpec,
        cds_spec::{CdsSpecMap, load_cds_spec},
        references::load_references,
    },
    data::{
        keys::{RefKey, SpecKey},
        weights::{CodonWeightMatrix, load_codon_weights},
    },
    toml::{AlignmentMethod, AlignmentParams, Formatting, Rules, TomlConfig},
};
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};
use zoe::{
    alignment::{Alignment, SharedProfiles},
    data::err::ResultWithErrorContext,
    iter_utils::ProcessResultsExt,
    prelude::*,
};

/// A fully materialized annotation module ready for protein annotation work.
#[derive(Debug)]
pub struct AnnotationModule<'a> {
    /// The name of the module (e.g., `flu`, `cov`, or `rsv`). This must
    /// correspond to a folder in `ribosome_res`.
    pub name:                    &'a String,
    /// An optional version for the module (e.g., `2.0-alpha`).
    pub version:                 &'a str,
    /// The method to use for performing sequence alignment.
    pub(crate) alignment_method: AlignmentMethod,
    pub(crate) formatting:       &'a Formatting,
    pub(crate) rules:            &'a Rules,
    /// A vector of the other module names and the paths to their reference
    /// files, used in for providing warning messages when unrecognized compound
    /// types are encountered.
    pub(crate) other_modules:    Vec<(&'a String, PathBuf)>,
    /// Compound type map for iteration-based processing.
    ///
    /// The keys of the map (ctypes) will not contain tabs.
    pub(crate) ctype_map:        HashMap<String, Vec<ReferenceGroup<'a>>>,
    /// Do we have codon-position weights to work with
    pub(crate) have_weights:     bool,
}

impl<'a> AnnotationModule<'a> {
    /// Creates a new [`AnnotationModule`] from the parsed TOML file by reading
    /// the references, specs, and other module files.
    ///
    /// ## Panics
    ///
    /// The `toml_path` must have a parent.
    ///
    /// ## Errors
    ///
    /// - The requested `module_name` must be present in the config.
    /// - The reference file in the config must exist within the module's folder
    ///   and be parsed successfully.
    /// - The codon position weights file if specified must exist within the
    ///   module's folder and be parsed successfully.
    /// - The CDS specifications file must exist within the module's folder and
    ///   be parsed successfully.
    pub fn new(config: &'a TomlConfig, toml_path: &Path, module_name: &str) -> std::io::Result<AnnotationModule<'a>> {
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

        // Load the reference sequences
        let references = {
            let references_path = module_root.join(&module.references);
            load_references(&references_path)
                .with_path_context("Failed to load the references from file", &references_path)?
        };

        // Load the codon weights if specified in TOML
        let mut codon_weights = if let Some(weights) = &module.weights {
            let weights_path = module_root.join(weights);
            load_codon_weights(&weights_path)
                .with_path_context("Failed to load the codon position weights from file", weights_path)?
        } else {
            CodonWeightMatrix::new()
        };

        let have_weights = !codon_weights.is_empty();

        // Load the CDS specs
        let mut cds_spec = {
            let cds_spec_path = module_root.join(&module.cds_spec);
            load_cds_spec(&cds_spec_path).with_path_context("Failed to load the CDS specs from file", cds_spec_path)?
        };

        let mut ctype_map: HashMap<String, Vec<ReferenceGroup>> = HashMap::new();

        for (ref_key, seqs) in references {
            let params = module.alignment.get(&ref_key.compound_type);

            // Get the list of groups for the given compound type
            let groups = ctype_map.entry(ref_key.compound_type.to_string()).or_default();

            // See if there is an existing entry in the list of groups for the given
            // reference ID
            if let Some(group) = groups.iter_mut().find(|group| group.reference_id == ref_key.reference_id) {
                // Update that reference group
                group.extend(seqs, &ref_key, params)?;
            } else {
                // Add a new reference group. Validity: seqs is non-empty by
                // guarantees on load_references
                groups.push(ReferenceGroup::new(
                    &ref_key,
                    seqs,
                    params,
                    &mut cds_spec,
                    &mut codon_weights,
                )?);
            }
        }

        Ok(AnnotationModule {
            ctype_map,
            name: &module.name,
            version: module.version.as_ref().map(AsRef::as_ref).unwrap_or_else(|| "unknown"),
            alignment_method: module.alignment_method,
            formatting: &module.formatting,
            rules: &module.rules,
            other_modules,
            have_weights,
        })
    }

    /// Do we have codon-position weights to work with
    pub fn have_weights(&self) -> bool {
        self.have_weights
    }

    /// Finds the best local Smith-Waterman alignment for a query sequence
    /// against all profiles in this group.
    ///
    /// The alignment with the highest score is considered best, or `None` is
    /// returned if no alignment was found. The `states` in the alignment will
    /// only include `M`, `I`, `D`, and `S`. Furthermore, the alignment is
    /// guaranteed to start and end with `M` states excluding soft clipping.
    pub(crate) fn best_alignment<T: AsRef<[u8]> + ?Sized>(
        &self, refs: &ReferenceGroup, query: &T,
    ) -> Option<Alignment<u32>> {
        // TODO: The use of get feels questionable. Even if it is
        // unlikely to ever happen, it feels like it should be an error/warning.
        let mut alignments = refs.profiles.iter().filter_map(|p| {
            let alignment = match self.alignment_method {
                // Validity: only contains M, I, D, and S per TODO
                AlignmentMethod::OnePass => p.sw_align_from_i16(query.as_query_src()),
                // Validity: only contains M, I, D, and S per TODO
                AlignmentMethod::ThreePass => p.sw_align_from_i16_3pass(query.as_query_src()),
            };

            alignment.get()
        });

        let mut best_alignment = alignments.next()?;

        for alignment in alignments {
            // Instead of using max_by_key, we use manual comparison to ensure
            // the first maximum is returned (not the last)
            if alignment > best_alignment {
                best_alignment = alignment;
            }
        }

        // Validity: the alignment cannot begin or end with indels since
        // ReferenceGroup guarantees that either gap_open or gap_extend are
        // strictly negative
        Some(best_alignment)
    }

    /// Attempts to return the module name of a different module containing the
    /// specified `ctype`.
    pub fn find_in_other_module(&self, ctype: &str) -> Option<&String> {
        for (module_name, ref_path) in &self.other_modules {
            if let Ok(reader) = FastaReader::from_path(ref_path) {
                for fasta_result in reader.flatten() {
                    if let Ok(key) = RefKey::parse(&fasta_result.name)
                        && ctype == key.compound_type
                    {
                        return Some(module_name);
                    }
                }
            }
        }

        None
    }
}

/// Pre-computed alignment profile for a reference sequence.
pub type AlignmentProfiles<'a> = SharedProfiles<'a, 32, 16, 8, 5>;

/// Information about references sharing the same `reference_id` within a
/// compound type.
///
/// All the references must be the same length.
#[derive(Debug)]
pub(crate) struct ReferenceGroup<'a> {
    /// The shared reference ID of the reference sequences.
    ///
    /// This field will not contain tabs.
    pub(crate) reference_id:  String,
    /// The shared length of the reference sequences.
    pub(crate) length:        usize,
    /// The alignment profiles corresponding to the sequences.
    ///
    /// At least one of `gap_open` and `gap_extend` will be strictly negative
    /// for all profiles, to ensure that optimal local alignments do not begin
    /// or end with indels.
    pub(crate) profiles:      Vec<AlignmentProfiles<'a>>,
    /// The specifications for all the protein products that can be formed from
    /// the references.
    pub(crate) product_specs: Vec<ProductSpec>,
}

impl<'a> ReferenceGroup<'a> {
    /// Initializes a new [`ReferenceGroup`] containing profiles for the
    /// specified sequences.
    ///
    /// Entries are removed from `cds_spec` and `codon_weights` and used to
    /// create the [`ProductSpec`] entries held by the [`ReferenceGroup`].
    ///
    /// ## Panics
    ///
    /// Out-of-bounds indexing will occur if `seqs` is empty.
    ///
    /// ## Errors
    ///
    /// The lengths of the sequences must all be the same and be non-empty, and
    /// the profiles must build successfully.
    fn new(
        ref_key: &RefKey, seqs: Vec<Nucleotides>, params: &'a AlignmentParams, cds_spec: &mut CdsSpecMap,
        codon_weights: &mut CodonWeightMatrix,
    ) -> std::io::Result<Self> {
        let length = seqs[0].len();
        for seq in &seqs {
            if seq.len() != length {
                return Err(std::io::Error::other(format!(
                    "Inconsistent reference lengths for '{reference_id}|{compound_type}'",
                    reference_id = ref_key.reference_id,
                    compound_type = ref_key.compound_type
                )));
            }
        }

        // Validity: gap_open and gap_extend restrictions are met due to
        // AlignmentParams guarantees
        let profiles = seqs
            .into_iter()
            .map(|seq| {
                AlignmentProfiles::new(seq, &params.matrix, params.gap_open, params.gap_extend)
                    .with_context("Failed to build alignment profiles")
            })
            .collect::<Result<Vec<_>, _>>()?;

        let product_specs = cds_spec
            .remove(ref_key)
            .into_iter()
            .flatten()
            .map(|(product_name, exons)| {
                let spec_key = SpecKey::new(&ref_key.reference_id, &product_name);
                ProductSpec {
                    name: product_name,
                    exons,
                    codon_weights: codon_weights.remove(&spec_key),
                }
            })
            .collect();

        Ok(Self {
            reference_id: ref_key.reference_id.to_string(),
            length,
            profiles,
            product_specs,
        })
    }

    /// Adds additional sequences to the [`ReferenceGroup`].
    ///
    /// ## Errors
    ///
    /// The lengths of the sequences must equal the lengths of the existing
    /// sequences in the [`ReferenceGroup`], and the profiles must build
    /// successfully.
    fn extend(&mut self, seqs: Vec<Nucleotides>, ref_key: &RefKey, params: &'a AlignmentParams) -> std::io::Result<()> {
        for seq in &seqs {
            if seq.len() != self.length {
                return Err(std::io::Error::other(format!(
                    "Inconsistent reference lengths for '{reference_id}|{compound_type}'",
                    reference_id = ref_key.reference_id,
                    compound_type = ref_key.compound_type
                )));
            }
        }

        // Validity: gap_open and gap_extend restrictions are met due to
        // AlignmentParams guarantees
        let profiles = seqs.into_iter().map(|seq| {
            AlignmentProfiles::new(seq, &params.matrix, params.gap_open, params.gap_extend)
                .with_context("Failed to build alignment profiles")
        });
        profiles.process_results(|iter| self.profiles.extend(iter))?;

        Ok(())
    }
}
