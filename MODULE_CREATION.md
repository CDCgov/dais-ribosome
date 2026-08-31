# Module Creation

A module in DAIS-ribosome, used to model different viruses, is represented by a
collection of files and an entry in a TOML file. The following steps can be used
to add a new module to DAIS-ribosome:

1. Choose a name for the module. This is advised to be a short, lowercase
   description, such as `flu`, `cov`, or `rsv`.
2. Create a folder with that name under `ribosome_res/`.
3. Determine the coordinate spaces, types, and proteins the module will model
   (see [Planning Your Module](#planning-your-module)).
4. Within the new folder, add three files:
    - A FASTA file with the reference sequences, as described in
      [Reference Sequence File](#reference-sequence-file).
    - A TSV file with the coding sequence specifications, as described in
      [Coding Sequence Specifications](#coding-sequence-specifications).
    - A TSV file with the codon-position weights, as described in
      [Codon-Position Weights](#codon-position-weights).
5. Add the module to `ribosome_res/modules.toml` according to
   [Module TOML Format](#module-toml-format).

Unlike `v1` which required a separate "build" step and the use of some perl
scripts, `v2` reads the files directly and does not require any other steps
besides those listed above.

## Planning Your Module

Before creating the files in the following steps, it is useful to plan out the
structure of the module.

A DAIS-ribosome module can model multiple related classifications, as well as
the different segments of a segmented genome. This information is combined into
a _compound type_ or _ctype_, which is any classification level so long as it is
internally consistent. If the genome is segmented, this must represent the
segment as well. For `flu`, the ctype combines the type, subtype, and segment,
such as `A_HA_H7`. For `cov`, which does not have a segmented genome, this is
currently either `SARS-CoV-2` or `MERS-CoV`.

To each ctype, there should be associated one or more references, grouped under
one or more _reference IDs_. In `cov` and `rsv`, only a single reference
sequence is used for each `ctype`, and this is typically a good starting place.
For a segmented genome, the reference sequences for each segment should be
grouped under the same reference ID.

If a single reference sequence is insufficient to model all the variation in a
given `ctype`, additional reference sequences can be added under each reference
ID. When aligning a sequence with a given ctype against a reference ID, all
pairwise alignments are performed and the highest scoring is used. Note that the
references grouped together under an ID should represent the same type, and for
a given segment, all the references within the ID must model the same coordinate
space and have the same length.

Additionally, adding additional reference IDs can allow for more outputs to be
formed for a single `ctype`. For example, in `flu`, B/Victoria and B/Yamagata
have relatively high homology and are not represented as different ctypes.
Instead, two reference IDs are included (`BRISBANE60` for Victoria and
`PHUKET3073` for Yamagata), and every `B_HA` and `B_NS` input produces outputs
as aligned against `BRISBANE60` _and_ as aligned against `PHUKET3073`.
Downstream programs may then perform classification or decide which records to
use.

Each `ctype` may produce multiple _protein products_, which can represent genes,
peptides, subunits of a protein, or any other coding region of interest. For
example, the `A_HA_H1` _ctype_ in `flu` produces `HA-signal`, `HA`, and `HA1`
products but the stalk protein `HA2` was purposefully omitted. Moreover, in the
`cov` module the `orf1ab` product is specified instead of the 16 non-structural
proteins that `orf1ab` encodes after processing.

Finally, since a product belongs to a single segment, and a reference ID
corresponds to a single taxonomic set of types, every pair of product and
reference ID should correspond to a single `ctype`. For example, reference
`CALI07` is from an influenza A(H1N1)pdm09 virus and product `HA1` is a subunit
of segment 4 (HA), so together `CALI07` and `HA1` map directly to `A_HA_H1`.

## Reference Sequence File

The reference sequences for a module are stored in a FASTA file. The headers
should be pipe-delimited in the form `<reference_id>|<compound_type>`. Any
additional pipe-delimited fields are ignored, but can be used to store other
metadata.

Both single- and multi-line FASTA files are supported.

During parsing, the sequences are validated to ensure they solely contain IUPAC
bases without gaps.

Below is an example. For a full example, see the [flu
module](ribosome_res/flu/flu-references.fasta).

```fasta
>CALI07|A_HA_H1|FJ981613
ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTGCAACCGCAAA...
>CALI07|A_HA_H1|FJ981613
ATGAAAGCAAAACTACTAGTCCTGTTATGTGCATTTACAGCTACATA...
>CALI07|A_NA_N1|GQ377078
ATGAATCCAAACCAAAAGATAATAACCATTGGTTCGGTCTGTATGAC...
>CALI07|A_NA_N1|GQ475920_DEL
ATGAACCCAAATCAAAAGATAATAACCATTGGATCAATCAGTATAGC...
>HK4801|A_PA|PZ114221
ATGGAAGATTTTGTGCGACAATGCTTCAACCCGATGATTGTCGAACT...
>HK4801|A_HA_H3|PZ114225
ATGAAGACTATCATTGCTTTGAGCTACATTCTATGTCTGGTTTTCGC...
>HK4801|A_PB1|PZ114223
ATGGATGTCAATCCGACTCTACTGTTCTTAAAAGTTCCAGCGCAAAA...
...
```

## Coding Sequence Specifications

The coding sequences (and their protein products) are specified in a TSV file
with the following columns:

1. `Reference ID`: The reference ID, which must match the reference ID specified
   in the reference sequence file.
2. `Compound Type`: The compound type (or ctype), which must match the compound
   type specified in the reference sequence file. This can be any classification
   level so long as it is internally consistent.
3. `Protein`: The name of the protein product.
4. `Coords`: The coordinates of the coding sequence, consisting of one or more
   ranges. The syntax for a single range is `<start>..<end>`, where both are
   1-based and inclusive. Multiple ranges are separated by semicolons. See the
   following [rules](#rules-for-exoncds-ranges) for more details.
5. `Required Beginning`: An optional field for a required first codon that must
   be present at the start of the product in order for it to be included. If
   present, this is often `ATG`. This can be useful for modeling products which
   may or may not be produced by a given sequence, depending on whether
   particular mutations are present. `PB1-F2` is an example for influenza.

Below is an example of the TSV contents. For a full example, see the [flu
module](ribosome_res/flu/flu-cds-spec.tsv).

| reference_id | ctype   | product_name | coords          | required_beginning |
| ------------ | ------- | ------------ | --------------- | ------------------ |
| CALI07       | A_HA_H1 | HA-signal    | 1..51           |                    |
| CALI07       | A_HA_H1 | HA           | 52..1701        |                    |
| CALI07       | A_HA_H1 | HA1          | 52..1032        |                    |
| CALI07       | A_NA_N1 | NA           | 1..1410         |                    |
| HK4801       | A_PA    | PA           | 1..2151         |                    |
| HK4801       | A_PA    | PA-X         | 1..570;572..760 |                    |
| HK4801       | A_PB1   | PB1          | 1..2274         |                    |
| HK4801       | A_PB1   | PB1-F2       | 95..2272        | ATG                |
| HK4801       | A_HA_H3 | HA-signal    | 1..48           |                    |
| HK4801       | A_HA_H3 | HA           | 49..1701        |                    |
| HK4801       | A_HA_H3 | HA1          | 49..1035        |                    |
| ...          | ...     | ...          | ...             | ...                |

### Rules for Exon/CDS Ranges

The length of all the ranges in the coordinates field must sum to a multiple of
three (i.e., it must contain complete codons). If there are multiple coding
sequence ranges, then the ranges must be in order, and adjacent ranges must
either have a noncoding region between them or must overlap by no more than 2
bases. If one range starts immediately after another range, then the two ranges
should be combined.

A noncoding region between two ranges is useful for modeling introns,
alternative splicing, or positive ribosomal frameshifts. For influenza, `NEP`
and `M2` are examples of alternative splicing, and `PA-X` is an example of a
positive ribosomal frameshift.

An overlap between two adjacent ranges is useful for modeling negative ribosomal
frameshifts. In COVID, `orf1ab` is an example of this.

## Codon-Position Weights

The codon-position weights file is a multi-section TSV file containing a TSV
comment to denote the start of a new section. They are of the form
`#<reference_id>|<product_name>`. Any extra pipe-delimited fields are ignored.

Within each section, the codon-position weights are specified in a three column
TSV:

1. The codon number (or amino acid position) represented as a 1-based integer.
   Positions are with respect to the coordinate space for the reference and
   product listed in the section header.
2. The complete codon, containing solely IUPAC characters (and no gaps).
3. The observation count (or relative weight) of the codon, with higher values
   indicating a greater preference for that site.

The codons are converted to uppercase, and `U` is replaced with `T`. The
codon/position pairs within each section must be unique.

To save space, consider filtering codons with minimal observations or weights
(e.g., singletons).

Below is an example. For a full example, see the [flu
module](ribosome_res/flu/flu-codon-position-weights.tsv).

```tsv
#CALI07|HA
1	GAC	19768
1	GAT	9
2	ACA	19653
2	ACT	44
2	AAA	29
2	ACC	23
3	TTA	17880
3	ATA	1766
3	TTG	51
...

#CALI07|HA-signal
1	ATG	19768
2	AAG	17800
...

#HK4801|HA
1	CAA	29396
1	CAG	101
...
```

## Module TOML Format

The TOML file contains configuration and file paths for all the modules. Each
module has fields and sections in the format:

```toml
[[module]]
name = "module_name"
references = "module-references.fasta"
weights = "module-codon-position-weights.tsv"
cds_spec = "module-cds-spec.tsv"

[module.alignment]
default = { match = 14, mismatch = 1, gap_open = 40, gap_extend = 1 }
```

The `[module.formatting]` section can be included to further configure the
output formatting, and the `[module.rules]` section can be used to configure
rules used during materialization. The alignment weights can be configured by
ctype by adding additional lines under `[module.alignment]`:

```toml
ctype = { match = 10, mismatch = 1, gap_open = 50, gap_extend = 1 }
```

The `[module.experimental]` section can be included to utilize experimental
features. These features are in active development or are under evaluation, and
therefore may be changed or removed at any time.

The following is a list of the field paths that can be specified:

|                                                   Field Path | Type                       | Required | Default      | Description                                                                                                                             |
| -----------------------------------------------------------: | -------------------------- | -------- | ------------ | --------------------------------------------------------------------------------------------------------------------------------------- |
|                                                         name | String                     | Yes      | -            | The name of the module, matching the name of the folder in ribosome_res                                                                 |
|                                                      version | String                     | No       | None         | The version for the module                                                                                                              |
|                                            alternative_names | List of Strings            | No       | Empty        | Any alternative names that can be used to refer to the module                                                                           |
|                                                   references | File Name                  | Yes      | -            | The name of the FASTA file containing the references, which should be located in module's folder                                        |
|                                                      weights | File Name                  | Yes      | -            | The name of the codon-position weights file, which should be located in the module's folder                                             |
|                                                     cds_spec | File Name                  | Yes      | -            | The name of the CDS specs file, which should be located in the module's folder                                                          |
|                                             alignment_method | "one-pass" or "three-pass" | No       | "three-pass" | The alignment method to use, with one-pass being faster for small genomes and three-pass being more memory efficient for larger genomes |
|                                      formatting.right_pad_aa | Boolean                    | No       | True         | Whether to add right padding to the aligned amino acid sequence                                                                         |
|                                     formatting.right_pad_cds | Boolean                    | No       | True         | Whether to add right padding to the aligned coding sequence                                                                             |
|                                     formatting.right_pad_gen | Boolean                    | No       | True         | Whether to add right padding to the genome alignment                                                                                    |
|                                     rules.try_stop_extension | Boolean                    | No       | False        | Enables the stop extension rule                                                                                                         |
|                                          rules.chew_to_start | Boolean                    | No       | False        | Enables the chew to start rule                                                                                                          |
|                                   rules.repairable_end_limit | Integer                    | No       | 0            | If non-zero, enables that number of bases to be re-added to either end of the local alignment                                           |
|           alignment.default.match<br />alignment.CTYPE.match | Integer                    | Yes      | -            | The score for a match, with optional CTYPE-specific overrides                                                                           |
|     alignment.default.mismatch<br />alignment.CTYPE.mismatch | Integer                    | Yes      | -            | The score/penalty for a mismatch, automatically converted to be non-positive, with optional CTYPE-specific overrides                    |
|     alignment.default.gap_open<br />alignment.CTYPE.gap_open | Integer                    | Yes      | -            | The score/penalty for opening a gap, automatically converted to be non-positive, with optional CTYPE-specific overrides                 |
| alignment.default.gap_extend<br />alignment.CTYPE.gap_extend | Integer                    | Yes      | -            | The score/penalty for extending a gap, automatically converted to be non-positive, with optional CTYPE-specific overrides               |
|                                   experimental.rewrite_rules | File Name                  | No       | None         | The name of the rewrite rules file, which should be located in the module's folder (*warning: this is experimental and may be changed*) |

See
[*Zoe*](https://docs.rs/zoe/0.0.29/zoe/alignment/sw/index.html#affine-gap-penalties)
for more details on the affine gap penalty.

## Sequence Classification

To support unclassified input sequences, a corresponding
[SSWSort](https://github.com/CDCgov/sswsort) module should be created with the
same name as the DAIS-ribosome module. Otherwise, the DAIS-ribosome module will
only work for classified input sequences or when `--assume-default-ctype` is
passed.

When using SSWSort, DAIS-ribosome requires that the `sswsort_res` folder is
located in the same directory as the `ribosome_res` folder. Otherwise, the
SSWSort modules may fail to be located, and an error may be issued for
unclassified input sequences at runtime. Note that the `install.sh` script
automatically handles loading the `sswsort_res` folder from the GitHub
repository.

## Experimental: Rewriting Rules

To achieve desired alignments, it is important to pick good references that are
close to the sequences being annotated, and to properly adjust the
codon-position weights.

However, there are still scenarios where these steps may fail to achieve the
desired results. For example, if a deletion occurs in an area with repeats in
the reference, then the genome alignment might arbitrarily pick which repeat to
count as the deletion. This choice is deterministic, so it is typically not an
issue, but if one placement is desired over another, then this can be
problematic.

One *experimental* solution to this is deletion rewriting. The semantics of the
file format and rule may change even without a minor version release, so use
with caution. Start by creating an additional TOML file within the module,
typically named `<module>-rewrite-rules.toml`. Then add rewrite rules according
to the following format:

|                         Field Path | Type                      | Required | Description                                                                             |
| ---------------------------------: | ------------------------- | -------- | --------------------------------------------------------------------------------------- |
|        genome.deletions.rule.ctype | String                    | Yes      | The compound type for which the rule is applied                                         |
| genome.deletions.rule.reference_id | String                    | Yes      | The reference ID for which the rule is applied                                          |
|         genome.deletions.rule.from | String (`<start>..<end>`) | Yes      | The 1-based end-inclusive genome range where if a deletion occurs, it will be rewritten |
|           genome.deletions.rule.to | String (`<start>..<end>`) | Yes      | The 1-based end-inclusive genome range to move the deletion to                          |

After performing genome alignment, DAIS-ribosome will check whether a deletion
occurs at any of the `from` ranges. If it does, then it will attempt to move the
deletion to the `to` range. These ranges must be the same length and not be
equal. An example is:

```toml
[[genome.deletions]]
rule = { ctype = "B_HA", reference_id = "BRISBANE60", from = "535..543", to = "529..537" }
```

The following is a thorough list of conditions that must hold for rewriting to
succeed:

- The `from` range is a deletion with exactly the specified length (it cannot be
  longer than the `from` range)
- The `to` range must be included in the alignment
- Positions in the `to` range that are not also in `from` must be matches
- Any positions between the two ranges must also be matches
- At least one additional match state must be present adjacent to the `to` range
  to prevent pathological cases

The rules are then applied from left to right based on the `from` range. Any of
these behaviors are subject to change.
