# Flu Module Notes

## References & Flu Subtypes

The flu module currently supports influenza A and B. Reference ID represent a
specific coordinate space for a well-known reference virus for the given gene
segments. While we sometimes provide multiple reference ID for certain gene
segments, we also add multiple references per reference ID with lineage-specific
structural mutations to help with query coverage.

### Influenza B

| Reference ID | Lineage  | Segments                         |
| ------------ | -------- | -------------------------------- |
| PHUKET3073   | Yamagata | PB1, PB2, PA, HA, NP, NA, MP, NS |
| BRISBANE60   | Victoria | HA, NS                           |

### Influenza A

| Reference ID  | Subtype / Lineage | Segments                         |
| ------------- | ----------------- | -------------------------------- |
| CALI07        | H1N1              | HA, NA                           |
| ANNARBOR60    | H2                | HA                               |
| HK4801        | H3N2              | PB2, PB1, PA, HA, NP, NA, MP, NS |
| VT1203        | H5                | HA                               |
| ANHUI01       | H7N9              | HA, NA                           |
| BGD0994       | H9                | HA                               |
| ONTARIO6118   | N4                | NA                               |
| RU1526        | N5 (Eurasian)     | NA                               |
| ALASKA4733    | N5 (American)     | NA                               |
| SICHUAN26221  | N6                | NA                               |
| NL219         | N7                | NA                               |
| ASTRAKHAN3212 | N8                | NA                               |

### Protein Products

All eight flu segments are included as compound types by the flu module. They generate the following products:

| Segment Name | Products                                      | Additional Details                                                                                                                                                                                                                                                                            |
| ------------ | --------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `PB2`        | `PB2`                                         | `PB2-S1` is not currently included.                                                                                                                                                                                                                                                           |
| `PB1`        | `PB1`, `PB1-F2`                               | `PB1-F2` is produced via leaky scanning and is a virulence factor. `PB1-N40` is also made via leaky scanning but is not currently included.                                                                                                                                                   |
| `PA`         | `PA`, `PA-X` (influenza A only)               | `PA-X` is produced via a +1 ribosomal frameshift and plays a role in host shutoff. `PA-N182` and `PA-N155` (with truncated N-terminus) are not currently included.                                                                                                                            |
| `HA`         | `HA-signal`, `HA`, `HA1`                      | `HA-signal` is the signal peptide. `HA` is the full protein before cleavage, after signal peptide removal. `HA1` is the globular head post-cleavage, used in protein modeling; its epitopes are of interest; only present for H1, H3, and flu B. `HA2` (the stalk) is not currently included. |
| `NP`         | `NP`                                          |                                                                                                                                                                                                                                                                                               |
| `NA`         | `NA`, `NB` (influenza B only)                 | `NB` is a flu B-specific product encoded by the NA segment via leaky scanning.                                                                                                                                                                                                                |
| `MP`         | `M1`, `M2` (influenza A), `BM2` (influenza B) | `M2` is alternatively spliced. `BM2` is produced via a stop-start pentanucleotide. The rare alternative splicing `M42` for influenza A is not included.                                                                                                                                       |
| `NS`         | `NS1`, `NEP`                                  | `NEP` is the alternatively-spliced nuclear export protein.                                                                                                                                                                                                                                    |

## Stop Codons

In previous releases of DAIS-ribosome, there was inconsistency with whether the `HA` products should include the stop codon or not. The flu module now includes stop codons for consistency, and to enable stop extension.
