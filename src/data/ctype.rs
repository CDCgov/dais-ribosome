//! Compound type hierarchy: ctype → reference_id → protein_product.
use crate::{
    config::{AlignmentParams, AlignmentWeights, product_spec::ProductSpec},
    data::{
        exons::Exons,
        keys::{RefKey, SpecKey},
        refs::ReferenceMap,
        spec::CdsSpecMap,
        weights::CodonWeightMatrix,
    },
};
use std::collections::HashMap;
use zoe::{
    alignment::{Alignment, SharedProfiles},
    data::{err::ResultWithErrorContext, nucleotides::Nucleotides},
    iter_utils::ProcessResultsExt,
    prelude::*,
};

/// Top-level index: ctype string → compound type data.
pub type CompoundTypeMap<'a> = HashMap<String, Vec<ReferenceGroup<'a>>>;

/// Pre-computed alignment profile for a reference sequence.
pub type AlignmentProfiles<'a> = SharedProfiles<'a, 32, 16, 8, 5>;

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
    pub fn new(
        ref_key: &RefKey, seqs: &'a [Nucleotides], params: &'a AlignmentParams,
        cds_spec: &mut HashMap<RefKey, Vec<(String, Exons)>>, codon_weights: &mut CodonWeightMatrix,
    ) -> std::io::Result<Self> {
        let length = seqs.first().map_or(0, Nucleotides::len);
        for seq in seqs {
            if seq.len() != length {
                return Err(std::io::Error::other(format!(
                    "Inconsistent reference lengths for '{reference_id}|{compound_type}'",
                    reference_id = ref_key.reference_id,
                    compound_type = ref_key.compound_type
                )));
            }
        }

        let profiles = seqs
            .iter()
            .map(|seq| {
                AlignmentProfiles::new(seq, &params.matrix, params.gap_open, params.gap_extend)
                    .with_context("Failed to build alignment profiles")
            })
            .collect::<Result<Vec<_>, _>>()?;

        let proteins = cds_spec
            .remove(ref_key)
            .unwrap_or_default()
            .into_iter()
            .map(|(protein_name, exons)| {
                let spec_key = SpecKey::new(&ref_key.reference_id, &protein_name);
                ProductSpec {
                    name: protein_name,
                    exons,
                    codon_weights: codon_weights.remove(&spec_key),
                }
            })
            .collect();

        Ok(Self {
            reference_id: ref_key.reference_id.to_string(),
            length,
            profiles,
            proteins,
        })
    }

    pub fn extend(&mut self, seqs: &'a [Nucleotides], ref_key: &RefKey, params: &'a AlignmentParams) -> std::io::Result<()> {
        for seq in seqs {
            if seq.len() != self.length {
                return Err(std::io::Error::other(format!(
                    "Inconsistent reference lengths for '{reference_id}|{compound_type}'",
                    reference_id = ref_key.reference_id,
                    compound_type = ref_key.compound_type
                )));
            }
        }

        let profiles = seqs.iter().map(|seq| {
            AlignmentProfiles::new(seq, &params.matrix, params.gap_open, params.gap_extend)
                .with_context("Failed to build alignment profiles")
        });
        profiles.process_results(|iter| self.profiles.extend(iter))?;

        Ok(())
    }

    /// Finds the best local Smith-Waterman alignment for a query sequence
    /// against all profiles in this group.
    ///
    /// The alignment with the highest score is considered best, or `None` is
    /// returned if no alignment was found. The `states` in the alignment will
    /// only include `M`, `I`, `D`, and `S`. Furthermore, any alignments are
    /// guaranteed to start and end with `M` states, excluding soft clipping.
    pub fn best_alignment<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<Alignment<u32>> {
        // TODO: The use of get feels questionable. Even if it is
        // unlikely to ever happen, it feels like it should be an error/warning.
        let mut alignments = self
            .profiles
            .iter()
            .filter_map(|p| p.sw_align_from_i16(query.as_query_src()).get());

        let mut best_alignment = alignments.next()?;

        for alignment in alignments {
            // Instead of using max_by_key, we use manual comparison to ensure
            // the first maximum is returned (not the last)
            if alignment > best_alignment {
                best_alignment = alignment;
            }
        }

        Some(best_alignment)
    }
}

/// Builds the compound type map from raw loaded data.
pub(crate) fn build_ctype_map<'a>(
    references: &'a ReferenceMap, mut cds_spec: CdsSpecMap, mut codon_weights: CodonWeightMatrix,
    alignment_weights: &'a AlignmentWeights,
) -> std::io::Result<CompoundTypeMap<'a>> {
    let mut ctype_map: HashMap<String, Vec<ReferenceGroup<'a>>> = HashMap::new();

    for (ref_key, seqs) in references {
        let params = alignment_weights.get(&ref_key.compound_type);

        // Get the list of groups for the given compound type
        let groups = ctype_map.entry(ref_key.compound_type.to_string()).or_default();

        // See if there is an existing entry in the list of groups for the given
        // reference ID
        if let Some(group) = groups.iter_mut().find(|group| group.reference_id == ref_key.reference_id) {
            // Update that reference group
            group.extend(seqs, ref_key, params)?;
        } else {
            // Add a new reference group
            groups.push(ReferenceGroup::new(ref_key, seqs, params, &mut cds_spec, &mut codon_weights)?);
        }
    }

    Ok(ctype_map)
}
