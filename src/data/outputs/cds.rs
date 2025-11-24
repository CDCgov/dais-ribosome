use crate::{
    config::Formatting,
    data::products::{ComputedDeletion, ComputedInsertion, ComputedProduct, Product},
};
use std::fmt::{self, Display};

/// A single row for `.seq` output.
pub struct SeqRow<'a> {
    pub id: &'a str,
    pub ctype: &'a str,
    pub ref_id: &'a str,
    pub(crate) product: &'a Product<'a>,
    pub(crate) computed_product: &'a ComputedProduct,
    pub(crate) formatting: &'a Formatting,
}

impl Display for SeqRow<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let leading = self.product.leading_cds_gap_len();
        let trailing = self.product.trailing_cds_gap_len();

        // Pre-allocate one padding string and slice it for all padding needs
        let pad = ".".repeat(leading.max(trailing));

        // TODO: Not sure if this will work for data with partial codons.
        let aa_lpad = &pad[..leading / 3];
        let aa_rpad = if self.formatting.right_pad_aa {
            &pad[..trailing / 3]
        } else {
            ""
        };
        let cds_lpad = &pad[..leading];
        let cds_rpad = if self.formatting.right_pad_cds { &pad[..trailing] } else { "" };

        let d = self.computed_product;

        write!(
            f,
            concat!(
                "{id}\t{ctype}\t{ref_id}\t{prot}\t{vh}",
                "\t{aa_seq}\t{aa_lpad}{aa_aln}{aa_rpad}",
                "\t{cds_id}\t{ins}\t{shift}",
                "\t{cds_seq}\t{c_lpad}{cds_aln}{c_rpad}",
                "\t{q_coords}\t{c_coords}"
            ),
            id = self.id,
            ctype = self.ctype,
            ref_id = self.ref_id,
            prot = self.product.product_spec.name,
            vh = d.variant_hash,
            aa_seq = d.aa_seq,
            aa_lpad = aa_lpad,
            aa_aln = d.aa_aln,
            aa_rpad = aa_rpad,
            cds_id = d.cds_id,
            ins = d.has_insertion,
            shift = d.has_shift_indel,
            cds_seq = d.cds_seq,
            c_lpad = cds_lpad,
            cds_aln = d.cds_aln,
            c_rpad = cds_rpad,
            q_coords = d.query_coords,
            c_coords = d.cds_coords,
        )
    }
}

/// A single row for `.ins` output.
pub struct InsRow<'a> {
    pub id:               &'a str,
    pub ctype:            &'a str,
    pub ref_id:           &'a str,
    pub protein:          &'a str,
    pub(crate) insertion: &'a ComputedInsertion,
}

impl Display for InsRow<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ins = self.insertion;
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.id,
            self.ctype,
            self.ref_id,
            self.protein,
            ins.upstream_aa,
            ins.inserted_nucleotides,
            ins.inserted_residues,
            ins.upstream_nt,
            ins.codon_shift,
        )
    }
}

/// A single row for `.del` output.
pub struct DelRow<'a> {
    pub id: &'a str,
    pub ctype: &'a str,
    pub ref_id: &'a str,
    pub protein: &'a str,
    pub(crate) computed_product: &'a ComputedProduct,
    pub(crate) deletion: &'a ComputedDeletion,
}

impl Display for DelRow<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let del = self.deletion;
        let d = self.computed_product;

        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.id,
            self.ctype,
            self.ref_id,
            self.protein,
            d.variant_hash,
            del.del_aa_start,
            del.del_aa_end,
            del.del_aa_len,
            del.in_frame,
            d.cds_id,
            del.del_cds_start,
            del.del_cds_end,
            del.del_cds_len,
        )
    }
}
