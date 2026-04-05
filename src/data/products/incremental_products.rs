use crate::{
    IteratorExt,
    data::{
        products::{ComputedDeletion, ComputedInsertion, Product},
        ranges::{CdsDeletionRange, CdsInsertionRange, CdsMatchRange, CdsStateRange, InsertionIdx},
    },
};
use std::{fmt::Display, ops::Range};
use zoe::prelude::{AminoAcids, Len, Nucleotides, Slice, Translate};

/// All fields of [`ComputedProduct`] that are populated incrementally via
/// iterating over the `product_ranges`.
///
/// None of the fields are populated with any ranges or residues past the first
/// encountered stop codon. From these fields, the rest of [`ComputedProduct`]
/// can be computed.
///
/// [`ComputedProduct`]: crate::data::products::ComputedProduct
pub struct ComputedIncrementalProducts {
    /// CDS alignment (with `-` for deletions, no insertions)
    pub cds_aln:                Nucleotides,
    /// Amino acid alignment (with `-` for deletions)
    pub aa_aln:                 AminoAcids,
    /// CDS sequence (without deletions, includes insertions)
    pub cds_seq:                Nucleotides,
    /// Whether any insertion exists in this product
    pub has_insertion:          bool,
    /// Whether any insertion or deletion causes a frameshift (length % 3 != 0)
    pub has_shift_indel:        bool,
    /// Query nucleotide coordinates (e.g., "1..45;48..90")
    pub query_coords:           Coords,
    /// CDS nucleotide coordinates (e.g., "1..45")
    pub cds_coords:             Coords,
    /// Computed non-filtered insertions for this product
    pub insertions:             Vec<ComputedInsertion>,
    /// Computed deletions for this product
    pub deletions:              Vec<ComputedDeletion>,
    /// The number of unaligned bases at the end of the coding sequence that
    /// were soft clipped or appeared after the first stop codon.
    ///
    /// This does not include trailing deletions, so that this field can be used
    /// to render right padding without double counting deletions.
    pub trailing_cds_unaligned: usize,
}

impl ComputedIncrementalProducts {
    /// Computes the incremental products from a `query` and `product`.
    pub fn new(query: &Nucleotides, product: &Product) -> Self {
        let range_capacity = product.product_ranges.len();

        let mut out = IncrementalAccumulator::new(product.leading_cds_unaligned, range_capacity);

        let end_cds_index = out.populate_from(query, product);

        let trailing_cds_unaligned = product.product_spec.exons.total_cds_length - end_cds_index;

        Self {
            cds_aln: out.cds_aln,
            aa_aln: out.aa_aln,
            cds_seq: out.dependent_fields.cds_seq,
            has_insertion: out.dependent_fields.has_insertion,
            has_shift_indel: out.dependent_fields.has_shift_indel,
            query_coords: out.dependent_fields.query_coords,
            cds_coords: out.dependent_fields.cds_coords,
            insertions: out.dependent_fields.insertions,
            deletions: out.dependent_fields.deletions,
            trailing_cds_unaligned,
        }
    }
}

/// A private struct for aiding in the construction of
/// [`ComputedIncrementalProducts`] which supports incremental updates for each
/// range.
struct IncrementalAccumulator {
    /// The aligned coding sequence (with `-` for deletions, no insertions).
    ///
    /// This is updated eagerly on each call to [`extend_from_match`] or
    /// [`extend_from_deletion`].
    ///
    /// [`extend_from_match`]: IncrementalAccumulator::extend_from_match
    /// [`extend_from_deletion`]: IncrementalAccumulator::extend_from_deletion
    cds_aln: Nucleotides,

    /// The aligned amino acid sequence, updated from `cds_aln` on each call to
    /// [`extend_from_match`] or when accumulation is complete. This is built
    /// incrementally to facilitate detection of the stop codon.
    ///
    /// [`extend_from_match`]: IncrementalAccumulator::extend_from_match
    aa_aln: AminoAcids,

    /// The index in `cds_aln` marking the start of the bases which have not yet
    /// been translated into `aa_aln`.
    ///
    /// If `untranslated_start == cds_aln.len()`, then all bases have been
    /// translated so far.
    untranslated_start: usize,

    /// All computed product fields which are dependent on the aligned
    /// sequences, and hence will be truncated if they are truncated.
    dependent_fields: DependentFields,
}

impl IncrementalAccumulator {
    /// Initializes an [`IncrementalAccumulator`] containing just the gap
    /// indicated by `leading_gap_len`.
    fn new(leading_cds_unaligned: usize, range_capacity: usize) -> Self {
        Self {
            cds_aln:            Nucleotides::from(vec![b'.'; leading_cds_unaligned]),
            aa_aln:             AminoAcids::new(),
            untranslated_start: 0,
            dependent_fields:   DependentFields::new(range_capacity),
        }
    }

    /// Populates the [`IncrementalAccumulator`] from all ranges in `product`.
    fn populate_from(&mut self, query: &Nucleotides, product: &Product) -> usize {
        for state in &product.product_ranges {
            match state {
                CdsStateRange::M(m) => {
                    // Short circuit if a stop codon is found
                    if let Some(end_cds_index) = self.extend_from_match(query, m) {
                        return end_cds_index;
                    }
                }
                CdsStateRange::I(ins) => self.extend_from_insertion(query, ins),
                CdsStateRange::D(del) => self.extend_from_deletion(del),
            };
        }

        // Add stop extension if the alignment includes one. This does not occur
        // if a stop codon was reached above.
        if let Some(ext_range) = &product.stop_extension_query_range {
            self.dependent_fields
                .extend_from_stop_extension(query, ext_range, product.product_spec.exons.total_cds_length);
        }

        // No stop codon reached, so finish the translation
        self.update_translation();

        // Possibly add a partial codon
        if self.untranslated_start < self.cds_aln.len() {
            self.aa_aln.push(b'~');
        }

        // Return the end_cds_index
        product
            .product_ranges
            .iter()
            .rev()
            .find_map(|s| match s {
                CdsStateRange::M(m) => Some(m.cds_range.end),
                CdsStateRange::D(d) => Some(d.cds_range.end),
                _ => None,
            })
            .unwrap_or(0)
    }

    /// Updates the accumulator with a match range. If a stop codon is reached,
    /// translation halts and the end index in the coding sequence is returned.
    #[must_use]
    fn extend_from_match(&mut self, query: &Nucleotides, range: &CdsMatchRange) -> Option<usize> {
        let slice = &query[range.query_range.clone()];
        self.cds_aln.extend_from_slice(slice);

        self.update_translation();

        // Check if we reached the stop codon
        if self.aa_aln.ends_with(b"*") {
            let range_until_stop = {
                let untranslated_len = self.cds_aln.len() - self.untranslated_start;
                let mut r = range.clone();
                r.cut_end(untranslated_len);
                r
            };

            // Shrink the aligned coding sequence
            self.cds_aln.shorten_to(self.untranslated_start);

            self.dependent_fields.extend_from_match(query, &range_until_stop);

            Some(range_until_stop.cds_range.end)
        } else {
            self.dependent_fields.extend_from_match(query, range);
            None
        }
    }

    /// Updates the accumulator with an insertion range.
    fn extend_from_insertion(&mut self, query: &Nucleotides, range: &CdsInsertionRange) {
        self.dependent_fields.extend_from_insertion(query, range);
    }

    /// Updates the accumulator with a deletion range.
    fn extend_from_deletion(&mut self, range: &CdsDeletionRange) {
        self.cds_aln.extend(std::iter::repeat_n(b'-', range.len()));
        self.dependent_fields.extend_from_deletion(range);
    }

    /// Updates the `aa_aln` translation to include as many of the bytes in
    /// `cds_aln` as possible. Only complete (length-3) codons are translated.
    fn update_translation(&mut self) {
        let untranslated = self.cds_aln.slice(self.untranslated_start..);
        self.aa_aln
            .extend(untranslated.to_aa_iter_exact().take_until_inclusive(|aa| *aa == b'*'));
        self.untranslated_start = self.aa_aln.len() * 3;
    }
}

/// All incrementally computed product fields which are dependent on the aligned
/// sequences, and hence will be truncated if they are truncated.
struct DependentFields {
    /// The unaligned coding sequence (without deletions, includes insertions).
    cds_seq: Nucleotides,

    /// Whether any insertion exists in this product
    has_insertion: bool,

    /// Whether any insertion or deletion causes a frameshift (length % 3 != 0)
    has_shift_indel: bool,

    /// Query nucleotide coordinates (e.g., "1..45;48..90")
    query_coords: Coords,

    /// CDS nucleotide coordinates (e.g., "1..45")
    cds_coords: Coords,

    /// Computed non-filtered insertions for this product
    insertions: Vec<ComputedInsertion>,

    /// Computed deletions for this product
    deletions: Vec<ComputedDeletion>,
}

impl DependentFields {
    /// Initializes an empty [`DependentFields`] with the given number of ranges
    /// as a starting capacity for [`Coords`] fields.
    fn new(range_capacity: usize) -> Self {
        Self {
            cds_seq:         Nucleotides::new(),
            has_insertion:   false,
            has_shift_indel: false,
            query_coords:    Coords::with_capacity((5 + 2 + 5) * range_capacity),
            cds_coords:      Coords::with_capacity((5 + 2 + 5) * range_capacity),
            insertions:      Vec::new(),
            deletions:       Vec::new(),
        }
    }

    /// Extends the dependent fields from a possibly-truncated match range.
    fn extend_from_match(&mut self, query: &Nucleotides, range: &CdsMatchRange) {
        let slice = &query[range.query_range.clone()];
        self.cds_seq.extend_from_slice(slice);

        self.query_coords.push_range(&range.query_range);
        self.cds_coords.push_range(&range.cds_range);
    }

    /// Extends the dependent fields from a insertion range.
    fn extend_from_insertion(&mut self, query: &Nucleotides, range: &CdsInsertionRange) {
        let slice = &query[range.query_range.clone()];

        // 0-based index after is equivalent to 1-based index
        // before. Validity: slice is from query, which meets
        // requirements
        let (computed_insertion, filtered) = ComputedInsertion::new(range.cds_index, slice);

        self.cds_seq.extend_from_slice(&computed_insertion.inserted_nucleotides);

        if !filtered {
            self.insertions.push(computed_insertion);
        }
        self.has_insertion = true;

        self.query_coords.push_range(&range.query_range);

        // If the insertion happens at the beginning of the sequence, do not
        // include this coordinate in cds_coords.
        if !range.cds_index.at_start() {
            self.cds_coords.push_upstream(range.cds_index);
        }

        if !range.len().is_multiple_of(3) {
            self.has_shift_indel = true;
        }
    }

    /// Extends the dependent fields from a deletion range.
    fn extend_from_deletion(&mut self, range: &CdsDeletionRange) {
        if !range.len().is_multiple_of(3) {
            self.has_shift_indel = true;
        }

        let in_frame = range.cds_range.start.is_multiple_of(3) && range.len().is_multiple_of(3);

        // TODO: this behavior will need regression tested

        // The 1-based inclusive start. Floor divide so that
        // deleting any part of a codon deletes the amino acid. Add
        // 1 to make it 1-based.
        let del_aa_start = (range.cds_range.start / 3) + 1;

        // The 1-based inclusive end. Ceiling divide so that
        // deleting any part of a codon deletes the amino acid.
        // 0-based exclusive end is equivalent to 1-based inclusive
        // end.
        let del_aa_end = range.cds_range.end.div_ceil(3);

        // TODO: Modified due to regression tests
        let del_aa_len = range.cds_range.len().div_ceil(3);

        self.deletions.push(ComputedDeletion {
            del_aa_start,
            del_aa_end,
            del_aa_len,
            in_frame,
            del_cds_start: range.cds_range.start + 1,
            del_cds_end: range.cds_range.end,
            del_cds_len: range.len(),
        });
    }

    /// Extends the dependent fields from a stop extension.
    fn extend_from_stop_extension(&mut self, query: &Nucleotides, ext_range: &Range<usize>, total_cds_length: usize) {
        if let Some(slice) = query.get(ext_range.clone()) {
            let nt_insertion_idx = InsertionIdx::from_right_idx(total_cds_length);
            // Validity: slice is from query, which meets validity requirements
            let (ins, filtered) = ComputedInsertion::new(nt_insertion_idx, slice);
            if !filtered {
                self.cds_seq.extend_from_slice(&ins.inserted_nucleotides);
                self.insertions.push(ins);
            }
            // TODO: Added based on regression tests:
            self.has_insertion = true;

            self.query_coords.push_range(ext_range);

            self.cds_coords.push_upstream(nt_insertion_idx);

            if !ext_range.len().is_multiple_of(3) {
                self.has_shift_indel = true;
            }
        }
    }
}

/// A helper struct for incrementally combining ranges/indices into a
/// [`String`].
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Coords(String);

impl Coords {
    /// Creates a new [`Coords`] such that the underlying [`String`] has the
    /// specified `capacity`.
    pub fn with_capacity(capacity: usize) -> Self {
        Self(String::with_capacity(capacity))
    }

    /// Returns the coordinates as a string reference.
    pub fn as_str(&self) -> &str {
        &self.0
    }

    /// PUshes a range to the coordinates.
    pub fn push_range(&mut self, range: &Range<usize>) {
        if !self.0.is_empty() {
            self.0.push(';');
        }

        let mut buff = core::fmt::NumBuffer::new();

        // 0-based half-open to 1-based inclusive
        self.0.push_str((range.start + 1).format_into(&mut buff));
        self.0.push_str("..");
        self.0.push_str(range.end.format_into(&mut buff));
    }

    // TODO: Rename?

    /// Pushes an insertion index to the coordinates.
    pub fn push_upstream(&mut self, index: InsertionIdx) {
        if !self.0.is_empty() {
            self.0.push(';');
        }

        let mut buff = core::fmt::NumBuffer::new();

        self.0.push_str(index.left_pos().format_into(&mut buff));
    }
}

impl Display for Coords {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl AsRef<str> for Coords {
    fn as_ref(&self) -> &str {
        &self.0
    }
}
