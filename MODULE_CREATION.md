# Module Creation

A module in DAIS-ribosome, used to model different viruses, is represented by a
collection of files and an entry in a TOML file. The following steps can be used
to add a new module to DAIS-ribosome:

1. Choose a name for the module. This is advised to be a short, lowercase
   description, such as `flu`, `cov`, or `rsv`.
2. Create a folder with that name under `ribosome_res/`.
3. Within the new folder, add three files:
    - A FASTA file with the reference sequences, as described in
      [Reference Sequence File](#reference-sequence-file).
    - A TSV file with the coding sequence specifications, as described in
      [Coding Sequence Specifications](#coding-sequence-specifications).
    - A TSV file with the codon-position weights, as described in
      [Codon-Position Weights](#codon-position-weights).
4. Add the module to `ribosome_res/modules.toml` according to
   [Module TOML Format](#module-toml-format).

Unlike `v1` which required a separate "build" step and the use of some perl
scripts, `v2` reads the files directly and does not require any other steps
besides those listed above.

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

```
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
|--------------|---------|--------------|-----------------|--------------------|
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

```
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

See
[*Zoe*](https://docs.rs/zoe/0.0.29/zoe/alignment/sw/index.html#affine-gap-penalties)
for more details on the affine gap penalty.
