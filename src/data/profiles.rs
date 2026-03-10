//! Alignment profiles and weight matrices.

use crate::{
    config::{AlignmentParams, ConfiguredModule},
    data::error::ModuleLoadError,
};
use std::collections::HashMap;
use zoe::{
    alignment::{ProfileError, SharedProfiles},
    data::{matrices::WeightMatrix, nucleotides::Nucleotides},
};

const DNA_ALPHA_LEN: usize = 5;

/// Collection of alignment weights for a module.
///
/// Provides default scoring parameters plus optional overrides for specific
/// compound types.
#[derive(Debug)]
pub struct AlignmentWeights {
    /// Default scoring matrix and gap parameters.
    pub default:   (AlignmentParams, WeightMatrix<'static, i8, DNA_ALPHA_LEN>),
    /// Per-compound-type overrides.
    pub overrides: HashMap<String, (AlignmentParams, WeightMatrix<'static, i8, DNA_ALPHA_LEN>)>,
}

impl AlignmentWeights {
    /// Build alignment weights from a module configuration.
    ///
    /// ## Errors
    ///
    /// The module must specify default alignment parameters.
    pub fn from_config(module: &ConfiguredModule) -> Result<Self, ModuleLoadError> {
        let default = module
            .alignment
            .get("default")
            .map(|params| {
                (
                    params.clone(),
                    WeightMatrix::new_dna_matrix(params.match_score, params.mismatch, Some(b'N')),
                )
            })
            .ok_or(ModuleLoadError::invalid_config(
                &module.name,
                "No default alignment parameters were specified",
            ))?;

        let mut overrides = HashMap::new();
        for (key, params) in &module.alignment {
            if key != "default" {
                let weights = (
                    params.clone(),
                    WeightMatrix::new_dna_matrix(params.match_score, params.mismatch, Some(b'N')),
                );
                overrides.insert(key.clone(), weights);
            }
        }

        Ok(Self { default, overrides })
    }

    /// Get the scoring parameters for a compound type.
    ///
    /// Returns the override for the compound type if it exists, otherwise returns defaults.
    pub fn get(&self, compound_type: &str) -> &(AlignmentParams, WeightMatrix<'static, i8, DNA_ALPHA_LEN>) {
        self.overrides.get(compound_type).unwrap_or(&self.default)
    }
}

/// Pre-computed alignment profile for a reference sequence.
pub type AlignmentProfiles<'a> = SharedProfiles<'a, 32, 16, 8, DNA_ALPHA_LEN>;

/// Build an alignment profile for a reference sequence.
pub fn build_alignment_profile<'a>(
    sequence: &'a Nucleotides, matrix: &'a WeightMatrix<'static, i8, DNA_ALPHA_LEN>, params: &AlignmentParams,
) -> Result<AlignmentProfiles<'a>, ProfileError> {
    SharedProfiles::new_with_w256(sequence.as_bytes(), matrix, params.gap_open, params.gap_extend)
}
