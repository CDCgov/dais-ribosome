//! Alignment profiles and weight matrices.

use crate::config::{AlignmentParams, ConfiguredModule};
use std::collections::HashMap;
use zoe::{
    alignment::SharedProfiles,
    data::{matrices::WeightMatrix, nucleotides::Nucleotides},
};

const DNA_ALPHABET: usize = 5;

/// Gap penalty parameters for alignment.
#[derive(Debug, Clone, Copy)]
pub struct GapParams {
    pub gap_open:   i8,
    pub gap_extend: i8,
}

impl GapParams {
    pub fn new(gap_open: i8, gap_extend: i8) -> Self {
        Self { gap_open, gap_extend }
    }
}

impl From<&AlignmentParams> for GapParams {
    fn from(params: &AlignmentParams) -> Self {
        Self {
            gap_open:   params.gap_open,
            gap_extend: params.gap_extend,
        }
    }
}

/// Collection of alignment weight matrices for a module.
///
/// Provides default scoring parameters plus optional per-compound-type overrides.
#[derive(Debug)]
pub struct AlignmentWeights {
    /// Default scoring matrix and gap parameters.
    pub default:   (GapParams, WeightMatrix<'static, i8, DNA_ALPHABET>),
    /// Per-compound-type overrides.
    pub overrides: HashMap<String, (GapParams, WeightMatrix<'static, i8, DNA_ALPHABET>)>,
}

impl AlignmentWeights {
    /// Build alignment weights from a module configuration.
    ///
    /// # Panics
    ///
    /// Panics if the module doesn't specify default alignment parameters.
    pub fn from_config(module: &ConfiguredModule) -> Self {
        let default = module
            .alignment
            .get("default")
            .map(|params| {
                (
                    GapParams::from(params),
                    WeightMatrix::new_dna_matrix(params.match_score, params.mismatch, Some(b'N')),
                )
            })
            .unwrap_or_else(|| panic!("Module '{}' must specify default alignment parameters", module.name));

        let mut overrides = HashMap::new();
        for (key, params) in &module.alignment {
            if key != "default" {
                let weights = (
                    GapParams::from(params),
                    WeightMatrix::new_dna_matrix(params.match_score, params.mismatch, Some(b'N')),
                );
                overrides.insert(key.clone(), weights);
            }
        }

        Self { default, overrides }
    }

    /// Get the scoring parameters for a compound type.
    ///
    /// Returns the override for the compound type if it exists, otherwise returns defaults.
    pub fn get(&self, compound_type: &str) -> &(GapParams, WeightMatrix<'static, i8, DNA_ALPHABET>) {
        self.overrides.get(compound_type).unwrap_or(&self.default)
    }
}

/// Pre-computed alignment profile for a reference sequence.
///
/// For reverse strand alignment, the query is reverse-complemented instead.
#[derive(Debug)]
pub struct AlignmentProfiles<'a> {
    pub profile: SharedProfiles<'a, 32, 16, 8, DNA_ALPHABET>,
}

impl<'a> AlignmentProfiles<'a> {
    /// Build a profile for a reference sequence.
    pub fn new(
        sequence: &'a Nucleotides, matrix: &'a WeightMatrix<'static, i8, DNA_ALPHABET>, params: &GapParams,
    ) -> Result<Self, super::error::ProfileBuildError> {
        let profile = SharedProfiles::new_with_w256(sequence.as_bytes(), matrix, params.gap_open, params.gap_extend)?;

        Ok(Self { profile })
    }
}
