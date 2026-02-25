use crate::{
    annotation::hashing::{nt_id, variant_hash},
    data::{
        products::{ComputedDeletion, ComputedInsertion, ComputedProduct, ProductSpec},
        ranges::{CdsMatchRanges, CdsStateRange},
    },
};
use std::{cmp::Ordering, iter::Extend, ops::Range};
use zoe::prelude::*;

#[derive(Debug)]
pub(crate) struct Product<'a> {
    /// The ranges within the exons that the query covers. This is initially
    /// formed by intersecting the query ranges with the exon ranges, then is
    /// tweaked.
    ///
    /// This is guaranteed to contain ordered and non-overlapping ranges. It
    /// does not begin or end with [`CdsStateRange::I`]. Within the aligned
    /// portion of the coding sequence (the exons), the ranges will be adjacent
    /// (forming a partition). However, there may be exons at the beginning or
    /// end which are not aligned against (or partially aligned against).
    pub(crate) product_ranges:             Vec<CdsStateRange>,
    pub(crate) product_spec:               &'a ProductSpec,
    /// If this product's last exon ends at the stop extension position, this
    /// holds the query range of the stop extension nucleotides.
    pub(crate) stop_extension_query_range: Option<Range<usize>>,
}

impl<'a> Product<'a> {
    /// Returns true only if the requirement exists and the codon did not match.
    pub(crate) fn missing_required_start(&self, query: impl AsRef<[u8]>) -> bool {
        if let Some(required) = self.product_spec.exons.required_start
            && let Some(CdsStateRange::M(m)) = self.product_ranges.first()
            && m.cds_range.start == 0
            && m.cds_range.end >= 3
            && let Some(mut first) = query
                .as_ref()
                .get(m.query_range.start..)
                .and_then(|s| s.first_chunk().copied())
        {
            first.make_ascii_uppercase();
            first != required
        } else {
            false
        }
    }

    pub(crate) fn leading_cds_gap_len(&self) -> usize {
        self.product_ranges
            .iter()
            .find_map(|s| match s {
                CdsStateRange::M(m) => Some(m.cds_range.start),
                CdsStateRange::D(d) => Some(d.cds_range.start),
                _ => None,
            })
            .unwrap_or(0)
    }

    /// Because the exon intersection algorithm can split deletions spanning
    /// multiple exons, this method merges adjacent deletion CDS for downstream
    /// correctness.
    pub(crate) fn condense_deletions(&mut self) {
        // The order of the closure arguments in dedup_by gets reversed.
        self.product_ranges.dedup_by(|range2, range1| {
            if let CdsStateRange::D(del1) = range1
                && let CdsStateRange::D(del2) = range2
            {
                // Ensure the deletions are adjacent in the coding sequence.
                // Validity: product_ranges forms an ordered partition over the
                // portion of the coding sequence which is aligned against.
                // Hence, adjacent deletions in product_ranges are guaranteed to
                // also be adjacent in the coding sequence.
                debug_assert_eq!(del1.cds_range.end, del2.cds_range.start);
                // Extend del1 (the kept element) to encompass del2
                del1.cds_range.end = del2.cds_range.end;
                // Return true to signal that deduplication is needed
                true
            } else {
                false
            }
        });
    }

    /// Fix indel frames by repositioning insertions and deletions to in-frame
    /// boundaries.
    ///
    /// This corrects in-frame indels (length divisible by 3) that occur at
    /// out-of-frame positions by shifting them left or right based on codon
    /// usage statistics. The algorithm mirrors the logic in
    /// `codonCorrectStats.pl`, albeit for coordinate math instead of strings.
    ///
    /// ## Validity
    ///
    /// `query_seq` should contain unaligned, uppercase IUPAC bases.
    pub(crate) fn fix_frames(&mut self, query: impl AsRef<[u8]>) {
        let len = self.product_ranges.len();
        let query_seq = query.as_ref();

        // We iterate by index because we need to access neighbors
        let mut i = 0;
        while i < len {
            match &self.product_ranges[i] {
                CdsStateRange::I(ins) => {
                    // Upstream index becomes upstream 1-based position
                    let frame = ins.frame();

                    // Only correct in-frame insertions at out-of-frame positions
                    if frame != 0 && ins.len() % 3 == 0 {
                        // Validity: this function requires the same
                        // requirements on query_seq
                        self.fix_insertion_frame(i, query_seq);
                    }
                }
                CdsStateRange::D(del) => {
                    // The inclusive 0-based index is the upstream 1-based position
                    let frame = del.cds_range.start % 3;

                    // Only correct in-frame deletions at out-of-frame positions
                    if frame != 0 && del.len() % 3 == 0 {
                        self.fix_deletion_frame(i, query_seq);
                    }
                }
                CdsStateRange::M(_) => {}
            }
            i += 1;
        }
    }

    // TODO: Can fixing frame mess up order of partition in case of flanking
    // deletion?

    /// Fixes an insertion at `idx` by shifting it to an in-frame position.
    ///
    /// Uses the A1/A2 shift logic from `codonCorrectStats.pl`:
    /// - A1 (frame 1): Compare codons formed by:
    ///   - 1 preceding base + first 2 insert bases
    ///   - last 1 insert base + 2 following
    ///+
    /// - A2 (frame 2): Compare codons formed by:
    ///   - 2 preceding bases + first 1 insert base
    ///   - last 2 insert bases + 1 following
    ///
    /// ## Validity
    ///
    /// `query_seq` should contain unaligned, uppercase IUPAC bases.
    fn fix_insertion_frame(&mut self, idx: usize, query_seq: &[u8]) -> Option<()> {
        let Product {
            product_ranges,
            product_spec,
            ..
        } = self;

        // Insertions cannot be split over exons, and in practice many scoring
        // schemes tend to avoid adjacent insertions and deletions, so in
        // practice left_match and right_match occur at idx-1 and idx+1
        // respectively. But in case either is a deletion, we perform a search.
        let Some((left_match, CdsStateRange::I(ins), right_match)) = Self::partition_states(product_ranges, idx) else {
            return None;
        };

        let frame = ins.frame();
        let codon_index = ins.codon_index();
        let insert_len = ins.len();

        // for weight lookup
        let codon_position = codon_index + 1;
        let insert_seq = &query_seq[ins.query_range.clone()];

        if frame == 2 {
            // A2: insertion after 2nd codon base
            //
            // A2/L2: shift insert left 2
            //   New codon: last 2 bases of insert + cp3
            //
            // A2/R1: shift insert right 1
            //   New codon: cp1 + cp2 + first 1 base of insert
            let cp1 = *query_seq.get(left_match.query_range.end - 2)?;
            let cp2 = *query_seq.get(left_match.query_range.end - 1)?;
            let cp3 = *query_seq.get(right_match.query_range.start)?;

            let a2l2 = [insert_seq[insert_len - 2], insert_seq[insert_len - 1], cp3];
            let a2r1 = [cp1, cp2, insert_seq[0]];

            // Validity: The codons are uppercase because they are derived from
            // bases in query_seq
            if product_spec.codon_left_ge_right(a2r1, a2l2, codon_position as u32) {
                // Insertion shifts right 1 for frame 2 split codon
                left_match.extend_end(1);
                ins.shift_right(1);
                right_match.shrink_start(1);
            } else {
                // Insertion shifts left 2 for frame 2 split codon
                left_match.shrink_end(2);
                ins.shift_left(2);
                right_match.extend_start(2);
            }
        } else if frame == 1 {
            // A1 insertion: insertion after 1st base of codon
            //
            // A1/L1: shift insert left 1
            //   New codon: last 1 base of insert + cp2 + cp3
            //
            // A1/R2: shift insert right 2
            //   New codon: cp1 + first 2 bases of insert

            let cp1 = *query_seq.get(left_match.query_range.end - 1)?;
            let cp2 = *query_seq.get(right_match.query_range.start)?;
            let cp3 = *query_seq.get(right_match.query_range.start + 1)?;

            let a1l1 = [insert_seq[insert_len - 1], cp2, cp3];
            let a1r2 = [cp1, insert_seq[0], insert_seq[1]];

            // Validity: The codons are uppercase because they are derived from
            // bases in query_seq
            if product_spec.codon_left_ge_right(a1l1, a1r2, codon_position as u32) {
                // Insertion shifts left 1 for frame 1 split codon
                left_match.shrink_end(1);
                ins.shift_left(1);
                right_match.extend_start(1);
            } else {
                // Insertion shifts right 2 for frame 1 split codon
                left_match.extend_end(2);
                ins.shift_right(2);
                right_match.shrink_start(2);
            }
        }

        Some(())
    }

    /// Fix a deletion at index `i` by shifting it to an in-frame position.
    fn fix_deletion_frame(&mut self, idx: usize, query: impl AsRef<[u8]>) -> Option<()> {
        let Product {
            product_ranges,
            product_spec,
            ..
        } = self;
        let query_seq = query.as_ref();

        let Some((left_match, CdsStateRange::D(del), right_match)) = Self::partition_states(product_ranges, idx) else {
            return None;
        };

        // Codon positions (1-based) at the boundaries of the deletion for table
        // checking
        let pos_left = (del.cds_range.start / 3) + 1;
        let pos_right = ((del.cds_range.end - 1) / 3) + 1;

        let frame = del.cds_range.start % 3;

        // Get the pivot codon bases that cross the deletion boundary
        if frame == 1 {
            // Frame 1 pivot codon: is last 1 base before gap + first 2 bases after gap
            let cp1 = query_seq.get(left_match.query_range.end - 1)?.to_ascii_uppercase();
            let cp2 = query_seq.get(right_match.query_range.start)?.to_ascii_uppercase();
            let cp3 = query_seq.get(right_match.query_range.start + 1)?.to_ascii_uppercase();

            let pivot = [cp1, cp2, cp3];

            // TODO: Likely want to compare with Ordering::is_gt instead, so
            // that ties resolve as right shift.

            // By default, we shift right for frame 1 (causes movement of 1
            // base, rather than 2). Only if there is evidence for left shift do
            // we shift left. Validity: The codon is uppercase because it is
            // derived from bases in query
            let shift_left = product_spec
                .compare_codon_positions(pos_left as u32, pos_right as u32, pivot)
                .is_some_and(Ordering::is_ge);

            if shift_left {
                // Deletion shift right 2 for frame 1 (2 bases move left)
                left_match.extend_end(2);
                del.shift_right(2);
                right_match.shrink_start(2);
            } else {
                // Deletion shift left 1 for frame 1 (1 base moves right)
                left_match.shrink_end(1);
                del.shift_left(1);
                right_match.extend_start(1);
            }
        } else if frame == 2 {
            // Frame2 pivot codon is last 2 bases before gap + first 1 base after gap.
            let cp1 = query_seq.get(left_match.query_range.end - 2)?.to_ascii_uppercase();
            let cp2 = query_seq.get(left_match.query_range.end - 1)?.to_ascii_uppercase();
            let cp3 = query_seq.get(right_match.query_range.start)?.to_ascii_uppercase();

            let pivot = [cp1, cp2, cp3];

            // By default, we shift the deletion right for frame 2 (causing 1
            // match to shift left, rather than 2). Only if there is evidence
            // for a left shift do we do that. pos_right being better means
            // shifting matches right is better, which means shifting deletion
            // left is better. Validity: The codon is uppercase because it is
            // derived from bases in query
            let shift_del_left = product_spec
                .compare_codon_positions(pos_left as u32, pos_right as u32, pivot)
                .is_some_and(Ordering::is_lt);

            if shift_del_left {
                // Deletion left shift 2 for frame 2 (2 bases move right)
                left_match.shrink_end(2);
                del.shift_left(2);
                right_match.extend_start(2);
            } else {
                // Deletion right shift 1 for frame 2 (1 base moves left)
                left_match.extend_end(1);
                del.shift_right(1);
                right_match.shrink_start(1);
            }
        }

        Some(())
    }

    /// Get the match ranges to modify and the indel
    fn partition_states(
        product_ranges: &mut [CdsStateRange], idx: usize,
    ) -> Option<(&mut CdsMatchRanges, &mut CdsStateRange, &mut CdsMatchRanges)> {
        let (left, mid, right) = product_ranges.split_around_mut(idx)?;

        let l = left.iter_mut().rev().filter_map(CdsStateRange::match_range_mut).next()?;
        let r = right.iter_mut().filter_map(CdsStateRange::match_range_mut).next()?;

        Some((l, mid, r))
    }

    pub(crate) fn add_query_coords(&mut self, offset: usize) {
        if offset == 0 {
            return;
        }

        for state in &mut self.product_ranges {
            match state {
                CdsStateRange::M(m) => m.shift_query_right(offset),
                CdsStateRange::I(i) => i.shift_query_right(offset),
                _ => {}
            }
        }
    }
}

/// An extension trait for slices, providing functionality specific to
/// DAIS-ribosome.
pub(crate) trait SliceExt<T> {
    /// Extracts a mutable index from a slice, along with the slice before it
    /// and the slice after it.
    fn split_around_mut(&mut self, index: usize) -> Option<(&mut [T], &mut T, &mut [T])>;
}

impl<T> SliceExt<T> for [T] {
    fn split_around_mut(&mut self, index: usize) -> Option<(&mut [T], &mut T, &mut [T])> {
        let (left, rest) = self.split_at_mut_checked(index)?;
        let (mid, right) = rest.split_first_mut()?;
        Some((left, mid, right))
    }
}

impl<'a> Product<'a> {
    /// Compute output data from this product.
    pub fn materialize(self, query: &Nucleotides) -> ComputedProduct<'a> {
        let query_bytes = query.as_bytes();
        let range_count = self.product_ranges.len();

        // Pre-pad CDS alignment for proper 5' reading frame
        let leading = self.leading_cds_gap_len();
        let mut cds_seq_bytes = Vec::new();
        let mut cds_aln_bytes = vec![b'.'; leading];
        let mut query_coords = String::with_capacity((5 + 2 + 5) * range_count);
        let mut cds_coords = String::with_capacity((5 + 2 + 5) * range_count);
        let mut insertions = Vec::new();
        let mut deletions = Vec::new();
        let mut has_shift_indel = false;

        let mut has_coords = false;
        for state in self.product_ranges.iter() {
            match state {
                CdsStateRange::M(m) => {
                    let slice = &query_bytes[m.query_range.clone()];
                    cds_seq_bytes.extend_from_slice(slice);
                    cds_aln_bytes.extend_from_slice(slice);

                    if has_coords {
                        query_coords.push(';');
                        cds_coords.push(';');
                    }
                    has_coords = true;
                    query_coords.push_range(&m.query_range);
                    cds_coords.push_range(&m.cds_range);
                }
                CdsStateRange::I(ins) => {
                    let slice = &query_bytes[ins.query_range.clone()];

                    cds_seq_bytes.extend_from_slice(slice);
                    if has_coords {
                        query_coords.push(';');
                        cds_coords.push(';');
                    }
                    has_coords = true;
                    query_coords.push_range(&ins.query_range);

                    // If the insertion happens at the beginning of the
                    // sequence, do not include this coordinate in cds_coords.
                    if ins.cds_index.index_after_ins() > 0 {
                        cds_coords.push_upstream(ins.cds_index.index_before_ins());
                    }

                    if !ins.len().is_multiple_of(3) {
                        has_shift_indel = true;
                    }

                    // 0-based index after is equivalent to 1-based
                    // index before
                    let upstream_nt = ins.cds_index.index_after_ins();
                    insertions.push(ComputedInsertion::new(upstream_nt, slice));
                }
                CdsStateRange::D(del) => {
                    cds_aln_bytes.extend(std::iter::repeat_n(b'-', del.len()));

                    let del_cds_len = del.len();
                    if !del_cds_len.is_multiple_of(3) {
                        has_shift_indel = true;
                    }

                    let in_frame = del.cds_range.start.is_multiple_of(3) && del_cds_len.is_multiple_of(3);

                    // TODO: this behavior will need regression tested
                    let del_aa_start = (del.cds_range.start / 3) + 1;
                    let del_aa_end = (del.cds_range.end - 1) / 3 + 1;
                    let del_aa_len = if in_frame {
                        del_cds_len / 3
                    } else {
                        del_aa_end - del_aa_start + 1
                    };

                    deletions.push(ComputedDeletion {
                        del_aa_start,
                        del_aa_end,
                        del_aa_len,
                        in_frame,
                        del_cds_start: del.cds_range.start + 1,
                        del_cds_end: del.cds_range.end,
                        del_cds_len,
                    });
                }
            }
        }

        // Compute has_insertion before adding stop extension inserts
        // so that the .seq output flag matches Perl makeProducts.pl behavior.
        let has_insertion = !insertions.is_empty();

        // If a stop extension applies to this product, materialize it as
        // a regular insertion at the end of the CDS.
        if let Some(ref ext_range) = self.stop_extension_query_range
            && let Some(slice) = query_bytes.get(ext_range.clone())
        {
            let upstream_nt = self.product_spec.exons.total_cds_length;
            let ins = ComputedInsertion::new(upstream_nt, slice);
            if !ins.filtered {
                insertions.push(ins);
            }
        }

        let cds_seq: Nucleotides = cds_seq_bytes.into();
        let cds_aln: Nucleotides = cds_aln_bytes.into();

        let aa_aln = cds_aln.translate_to_stop();

        // Derive aa_seq from aa_aln: strip alignment characters (`-` and `.`),
        // then splice in insertion residues at their correct AA positions.
        // This mirrors the Perl pipeline's `fa2delim -B -I` behavior and
        // ensures shift-deletions don't corrupt the reading frame.
        let mut aa_seq_bytes: Vec<u8> = aa_aln.iter().filter(|&&b| b != b'-' && b != b'.').copied().collect();

        let mut offset = 0;
        for ins in &insertions {
            let pos = ins.upstream_aa + offset;
            let residue_bytes = ins.inserted_residues.as_bytes();
            let insert_at = pos.min(aa_seq_bytes.len());
            for (j, &b) in residue_bytes.iter().enumerate() {
                aa_seq_bytes.insert(insert_at + j, b);
            }
            offset += residue_bytes.len();
        }

        let aa_seq = AminoAcids::from_vec_unchecked(aa_seq_bytes);
        let cds_id = nt_id(&cds_seq).unwrap_or_default();
        let variant_hash = variant_hash(&aa_seq).unwrap_or_default();

        ComputedProduct {
            product_ranges: self.product_ranges,
            product_spec: self.product_spec,
            cds_seq,
            cds_aln,
            cds_id,
            aa_seq,
            aa_aln,
            variant_hash,
            has_insertion,
            has_shift_indel,
            query_coords,
            cds_coords,
            insertions,
            deletions,
        }
    }
}

trait PushRange {
    fn push_range(&mut self, range: &Range<usize>);
    fn push_upstream(&mut self, index: usize);
}

impl PushRange for String {
    #[inline]
    fn push_range(&mut self, range: &Range<usize>) {
        let mut buff = core::fmt::NumBuffer::new();

        // 0-based half-open to 1-based inclusive
        self.push_str((range.start + 1).format_into(&mut buff));
        self.push_str("..");
        self.push_str(range.end.format_into(&mut buff));
    }

    #[inline]
    fn push_upstream(&mut self, index: usize) {
        let mut buff = core::fmt::NumBuffer::new();

        // 0-based upstream index to 1-based upstream position
        self.push_str((index + 1).format_into(&mut buff));
    }
}
