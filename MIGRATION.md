## Changes to Module Configuration/Specification Formats

In `v1`, the specifications for a module were represented by 4 files in `spec`:

- `<MODULE>.refs`, a FASTA file containing the references for the module
- `<MODULE>.spec`, the product specifications for the proteins
- `<MODULE>.ssw`, the weights to use when performing sequence alignment
- `<MODULE>.sto`, the codon/position weights in a perl storable file

In `v2`, the references, specifications, and codon/position weights are stored
in files under `ribosome_res/<MODULE>/`. The module is then defined in
`modules.toml`, where paths to the previous files and the sequence alignment
weights are specified. The TOML file also allows for customization of the rules
and formatting used by the module. 

Notably, the reference sequences no longer need to be reorganized into a `refs`
directory structure. This removes the need for the `rebuild` and makes editing
the reference sequences more direct.

In `v1`, the product specifications TSV file had rows of the form:

```
<Compound Type> \t <Protein> \t <Reference ID>|<Protein> \t <Coords> \t <Required Beginning>
```

Here, `|` represents a literal character and `<Required Beginning>` is optional.
`v2` changes this format to:

```
<Reference ID> \t <Compound Type> \t <Protein> \t <Coords> \t <Required Beginning>
```

The format of the coordinates field is now `start..end` where both are 1-based
inclusive. Previously, this was `start,end`.

The codon-position weights are now stored in a multi-section TSV file instead of
a perl storable file. Each section begins with a `#reference_id|protein` comment
line followed by records with position, codon, and count fields.

## Simplification of Built-in Modules

The modules have been renamed for conciseness. `INFLUENZA` is now `flu`,
`BETACORONAVIRUS` is now `cov`, and `RSV` is now `rsv`. However, the old module
names will continue to be recognized as alternative names.

The modules have also had all singleton codons (those with a count of 1) removed
from their codon-position weights to reduce the impact of noise (and the size of
the files).

## Padding Products, and Handling of Empty Products

When an alignment does not fully span all the exons of a product, padding (`.`)
is added to the start and/or end of the aligned sequence outputs. In the extreme
case where there is no intersection between the exons and the alignment, an
"empty product" is produced. Ribosome v2 has more consistent behavior for both
these cases:  

- In the case of an empty product, the previous version output `cds_aln` and
  `aa_aln` consisting of all `.`. The new version filters this output.
- The previous version never padded `aa_aln` on the right. Now it does so by
  default (customizable).
- If a product contains only a deletion, then `cds_seq` and other "empty" fields
  will contain the HADOOP null representation `\N` instead of being empty.

If a query sequence fails to align at all to a reference sequence (e.g., the
query is all `N`), then this output is filtered.

## Improvements to Frame Fixing (or Indel Shifting) Code

The frame fixing code is responsible for shifting out-of-frame indels to improve
the quality of the amino acid translation. Several improvements to the frame
fixing code were made to have better robustness and avoid edge cases:

- Deletions spanning or crossing into a non-coding region are no longer eligible
  to shift during frame fixing. Anecdotally this appears to enhance correctness.
- Indels within 3 bases of the end of a product are now eligible to shift.
- Any two deletions with less than 3 matched bases between them are now eligible
  to shift.
- Insertions at the end of a protein product or adjacent to a non-coding region
  are removed.

## Other Changes/Fixes from v1.7.1 to v2.0.0

- The CDS specifications for several HA products in the flu module have been
  extended to include the stop codon. This enables stop extension to occur, but
  will change may hashes/IDs.
- Nucleotides in the insertion output are no longer converted to lowercase.
- The warning "insertion exceeds range of annotated loci" is no longer produced.
- [*Zoe*]'s alignment implementation is used instead of `SSW`, so some
  alignments may change
- The query coordinates now have correct output. Certain patterns of indels in
  non-coding regions would previously cause incorrect query coordinate output.
- An insertion in the middle of a required start codon is now properly
  considered (it may cause the start codon to be interrupted, and hence the
  requirement not met).

[*Zoe*]: https://crates.io/crates/zoe