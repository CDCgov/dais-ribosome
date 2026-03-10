//! Alignment profiles.

use crate::config::toml::AlignmentParams;
use zoe::{
    alignment::{ProfileError, SharedProfiles},
    data::nucleotides::Nucleotides,
};

/// Pre-computed alignment profile for a reference sequence.
pub type AlignmentProfiles<'a> = SharedProfiles<'a, 32, 16, 8, 5>;

/// Build an alignment profile for a reference sequence.
pub fn build_alignment_profile<'a>(
    sequence: &'a Nucleotides, params: &'a AlignmentParams,
) -> Result<AlignmentProfiles<'a>, ProfileError> {
    SharedProfiles::new_with_w256(sequence.as_bytes(), &params.matrix, params.gap_open, params.gap_extend)
}
