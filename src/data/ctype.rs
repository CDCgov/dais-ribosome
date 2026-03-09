//! Compound type hierarchy: ctype → reference_id → protein_product.
use crate::{config::toml::AlignmentWeights, data::spec::CdsSpecMap};

use super::{
    error::ModuleLoadError,
    keys::{RefKey, SpecKey},
    products::ProductSpec,
    profiles::{AlignmentProfiles, build_alignment_profile},
    refs::ReferenceMap,
    weights::CodonWeightMatrix,
};
use std::collections::HashMap;
use zoe::alignment::Alignment;
use zoe::data::nucleotides::Nucleotides;
use zoe::prelude::*;

/// Top-level index: ctype string → compound type data.
pub type CompoundTypeMap<'a> = HashMap<String, Vec<ReferenceGroup<'a>>>;

/// Information about references sharing the same `reference_id` within a
/// compound type.
///
/// All the references must be the same length.
#[derive(Debug)]
pub(crate) struct ReferenceGroup<'a> {
    /// The shared reference ID of the reference sequences.
    pub(crate) reference_id: String,
    /// The shared length of the reference sequences.
    pub(crate) length:       usize,
    /// The alignment profiles corresponding to the sequences.
    pub(crate) profiles:     Vec<AlignmentProfiles<'a>>,
    pub(crate) proteins:     Vec<ProductSpec>,
}

impl<'a> ReferenceGroup<'a> {
    /// Finds the best alignment for a query sequence against all profiles in
    /// this group.
    ///
    /// Returns the alignment with the highest score, or `None` if no alignment
    /// was found.
    pub fn best_alignment<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<Alignment<u32>> {
        // TODO: Why did we hard code from i16? Will this potentially disagree
        // with aligner? Also, the use of get feels questionable. Even if it is
        // unlikely to ever happen, it feels like it should be an error/warning.
        self.profiles
            .iter()
            .filter_map(|p| p.sw_align_from_i16(query.as_query_src()).get())
            .max_by_key(|a| a.score)
    }
}

/// Build the ctype map from raw loaded data.
pub(crate) fn build_ctype_map<'a>(
    references: &'a ReferenceMap, mut cds_spec: CdsSpecMap, mut codon_weights: CodonWeightMatrix,
    alignment_weights: &'a AlignmentWeights,
) -> Result<CompoundTypeMap<'a>, ModuleLoadError> {
    // Regroup the sequences by compound_type and then reference_id (two levels
    // of grouping, rather than one using RefKey).
    let mut ctype_refs: HashMap<&str, Vec<(&str, &Vec<Nucleotides>)>> = HashMap::new();
    for (ref_key, seqs) in references {
        ctype_refs
            .entry(&ref_key.compound_type)
            .or_default()
            .push((&ref_key.reference_id, seqs));
    }

    let mut result = HashMap::new();

    for (ctype, ref_entries) in ctype_refs {
        let alignment_params = alignment_weights.get(ctype);

        // Group by reference_id within this ctype
        let mut ref_groups_map: HashMap<&str, Vec<&Nucleotides>> = HashMap::new();
        for (ref_id, seqs) in &ref_entries {
            ref_groups_map.entry(ref_id).or_default().extend(seqs.iter());
        }

        let mut reference_groups = Vec::with_capacity(ref_groups_map.len());

        for (reference_id, seqs) in ref_groups_map {
            // Ensure all sequences are the same length
            let length = seqs.first().map(|s| s.len()).unwrap_or(0);
            for seq in &seqs {
                if seq.len() != length {
                    return Err(ModuleLoadError::validation(format!(
                        "Inconsistent lengths for '{reference_id}|{ctype}'"
                    )));
                }
            }

            // Build the profiles for the equal-length sequences
            let profiles = seqs
                .iter()
                .map(|seq| build_alignment_profile(seq, alignment_params))
                .collect::<Result<Vec<_>, _>>()?;

            let proteins = cds_spec
                .remove(&RefKey::new(reference_id, ctype))
                .unwrap_or_default()
                .into_iter()
                .map(|(protein_name, exons)| {
                    let spec_key = SpecKey::new(reference_id, &protein_name);
                    ProductSpec {
                        name: protein_name,
                        exons,
                        codon_weights: codon_weights.remove(&spec_key),
                    }
                })
                .collect();

            reference_groups.push(ReferenceGroup {
                reference_id: reference_id.to_string(),
                length,
                profiles,
                proteins,
            });
        }

        result.insert(ctype.to_string(), reference_groups);
    }

    Ok(result)
}
