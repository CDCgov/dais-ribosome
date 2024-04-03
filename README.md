# DAIS-ribosome - annotation and coordinate mapping of CDS and protein sequences for virus genomes

The DAIS **ribosome** compartmentalizes the original translation engine developed for our protein analytics database. The tool was designed for use with INFLUENZA, but has been extended for use with BETACORONAVIRUS and RSV. I provide a brief outline of the algorithm:

1. If necessary, classify the nucleotide gene segment into its influenza type, segment, and subtype (compound or `C_type`) as well as make sure CDS is on the *forward strand*†.
2. Align said segments (via [SSW](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082138)) to the corresponding reference sequence(s) and pick the best alignment.
3. There can be more than one reference reading frame per C_type, so complete step 2 for each `reference_id`. Data with no reference is held-aside and added in later.
4. Fix alignment ends that have been chopped (due to local alignment disagreement).
5. Create protein product CDS using an internal specification.
6. Correct product alignment so that indels occur within frame only, then tabularize and create `cds_id`.
7. Amend insertion tables to use protein coordinates and translate. Filter singleton, doubleton, and entirely ambiugous inserts used for AA (retain in CDS output).
8. Translate CDS to amino acids, calculating the `variant_hash` as well.
9. Create coordinate mapping between CDS and AA
10. Combine AA, CDS, and coordinate tabular data; output with insertion data from step 6.
11. Produce deletion table(s) from aligned sequences.

**Table of contents:**

- [DAIS-ribosome - annotation and coordinate mapping of CDS and protein sequences for virus genomes](#dais-ribosome---annotation-and-coordinate-mapping-of-cds-and-protein-sequences-for-virus-genomes)
  - [Input data](#input-data)
  - [Output data](#output-data)
  - [Field descriptions](#field-descriptions)
  - [Special use of translated characters](#special-use-of-translated-characters)
  - [CLI Usage](#cli-usage)
  - [Installation](#installation)
  - [Updating the codon weight matrix](#updating-the-codon-weight-matrix)
  - [License](#license)
  - [Component Citations](#component-citations)

## Input data

**Input** can be one of the four formats.

1. Unannotated Fasta (ID only)

    >\>223550  
    >ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA

2. Annotated Fasta (ID and compound type)

    >\>223550|B_HA  
    >ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA

3. Unannotated tab-delimited (ID only)

    >223550 *\<tab\>* ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA

4. Annotated tab-delimited (ID and compound type)

    >223550 *\<tab\>* B_HA *\<tab\>* ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA

†Annotated fasta / tab-delimited data **must** be on the forward or plus strand. Unannotated data will be classified and reverse complemented to the forward strand in the CDS as needed.

## Output data

**Output** for the `dais-ribosome` consists of two tab-delimited files. One with `.ins` for insertions, one with `.del` for deletions, and `.seq` for sequence related data. An insertion file output example:

| ID | C_type | Ref_ID | Protein | Upstream_aa | Inserted_nuceotides | Inserted_residues | Upstream_nt | Codon_shift |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 11209 | B_HA | PHUKET3073 | HA | 161 | aaa | K | 483 | 0 |
| 154957 | B_HA | PHUKET3073 | HA | 163 | krc | X | 489 | 0 |
| 223550 | B_HA | PHUKET3073 | HA | 161 | caa | Q | 483 | 0 |

A deletion file example:

| ID | C_type | Ref_ID | Protein | VH | Del_AA_start | Del_AA_end | Del_AA_len | In_frame | CDS_ID | Del_CDS_start | Del_CDS_end | Del_CDS_len |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
|EPI_ISL_410721|SARS-CoV-2|WUHAN19|orf1ab|5ba70e95c9a3251bc6155f62295dd3e8|994|1002|9|true|29cd767e2d144c31179395fd606d1489ce731746|2980|3006|27|
|EPI_ISL_410721|SARS-CoV-2|WUHAN19|orf1ab|5ba70e95c9a3251bc6155f62295dd3e8|1012|1012|1|true|29cd767e2d144c31179395fd606d1489ce731746|3034|3036|3|
|EPI_ISL_410721|SARS-CoV-2|WUHAN19|S|450c068c437e7536d27fdb883d95d4f4|72|72|1|true|36a75a0d34960c048abaf82ee46a1b713eee534e|214|216|3|
|EPI_ISL_410721|SARS-CoV-2|WUHAN19|S|450c068c437e7536d27fdb883d95d4f4|146|146|1|true|36a75a0d34960c048abaf82ee46a1b713eee534e|436|438|3|
|EPI_ISL_410721|SARS-CoV-2|WUHAN19|S|450c068c437e7536d27fdb883d95d4f4|254|256|3|true|36a75a0d34960c048abaf82ee46a1b713eee534e|760|768|9|
|EPI_ISL_410721|SARS-CoV-2|WUHAN19|S|450c068c437e7536d27fdb883d95d4f4|680|683|4|true|36a75a0d34960c048abaf82ee46a1b713eee534e|2038|2049|12|

A sequence file output example:

| ID | C_type | Ref_ID | Protein | VH |  AA_seq | AA_aln | CDS_id | Insertion | Shift_Insert | CDS_seq | CDS_aln | Query_nt_coordinates | CDS_nt_coordinates |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 223550 | B_HA | BRISBANE60 | HA-signal | e81d2d895c70e91bb3ef917fe49fdab7 | MKAIIVLLMVVTSNA | MKAIIVLLMVVTSNA | 2aa6443b92ca45b301faa4d46e5fbd3b010e3ab7 |  false | false |ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | 20..64 | 1..45 |
| 223550 | B_HA | PHUKET3073 | HA-signal | e81d2d895c70e91bb3ef917fe49fdab7 | MKAIIVLLMVVTSNA | MKAIIVLLMVVTSNA | 2aa6443b92ca45b301faa4d46e5fbd3b010e3ab7 | false | false | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | 20..64 | 1..45 |
| 11209 | B_HA | BRISBANE60 | HA-signal | c7ee7ff234abf5c0591e0fe1af26ca87 | MKAIIILLMVVTSNA | MKAIIILLMVVTSNA | c49a73ab7280362c8c710abbf648708c41f97712 | false | false | ATGAAGGCAATAATTATACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTATACTACTCATGGTAGTAACATCCAATGCA | 1..45 | 1..45 |
| 11209 | B_HA | PHUKET3073 | HA-signal | c7ee7ff234abf5c0591e0fe1af26ca87 | MKAIIILLMVVTSNA | MKAIIILLMVVTSNA | c49a73ab7280362c8c710abbf648708c41f97712 | false | false | ATGAAGGCAATAATTATACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTATACTACTCATGGTAGTAACATCCAATGCA | 1..45 | 1..45 |

*Optional* genome output file example:

| ID | C_type | Ref_ID | Genome_ID | Genome_length | Insertion | Genome_seq | Genome_aln |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| EPI_ISL_FAKE1 | SARS-CoV-2 | WUHAN19 | 8e193a72b22a666947b21cb785af6780c2c6996b | 108 | true | TTTAAGGTTTATACCTTCCCAGGTAACAAACCAACC**TGGGTTTGG**AACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTC | .TTTAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTC |

*Optional* genome insertion output file example:

| ID | C_type | Ref_ID | Upstream_nt | Inserted_nucleotides |
| ------ | ------ | ------ | ------ | ------ |
| EPI_ISL_FAKE1 | SARS-CoV-2 | WUHAN19 | 37 | TGGGTTTGG |

*Optional* genome deletion output file example:

| ID | C_type | Ref_ID | Del_NT_start | Del_NT_end | Del_NT_len |
| ------ | ------ | ------ | ------ | ------ |  ------ |
| EPI_ISL_FAKE1 | SARS-CoV-2 | WUHAN19 |3246|3272|27|

## Field descriptions

The field explanations:

| field | description |
| ------ | ------ |
| ID | Any unique identifier, but likely the `flu_sequence_id`, `epi_segment_id`, raw `nt_id`, or NCBI accession. |
| C_type | The compound type conisting of the influenza type, segment, and subtype if applicable. This is the same as IRMA. Chimeric types start with an asterisk. Other modules, this field is used for the taxon, eg, SARS-CoV-2 and MERS-CoV. |
| Ref_ID | As with DAIS, the reference reading frame used for alignment. |
| Protein | The protein product or peptide derived from the gene segment. |
| VH (AA_ID) | The `variant_hash` as used in DAIS (md5 hex of `AA_seq`). |
| Insertion | Boolean indicating whether or not there is an insertion relative to reference in the original CDS. |
| Shift_Indel | Boolean indicating whether an insertion or deletion in the CDS is a non-triplet and could have induced a frameshift. |
| Upstream_aa / Upstream_nt | The upstream amino acid / nucleotide position for the insertion relative to the reference coordinates. |
| Codon | The codon(s) inserted. |
| Inserted_residues / Inserted_nucleotides | The residues / nucleotides inserted. |
| CDS_ID / Genome_ID | The nucleotide sequence ID using the sha1 hex of the CDS_seq and Genome_seq (as in PubSeq). |
| AA_seq / AA_aln | The amino acid sequence (less deletions + insertions) and the amino acid alignment (residues relative to reference only). |
| CDS_seq / CDS_aln | The nucleotide cds sequence (less deletions + insertions) and the cds alignment (bases relative to reference only). |
| Genome_seq / Genome_aln | The nucleotide genome sequence (less deletions + insertions) and the genome alignment (bases relative to reference only). |
| Query_nt_coordinates | Set of aligned position ranges representing the aligned coordinates relative to the original submitted query sequence. Insertions appear as singletons. |
| CDS_nt_coordinates | Set of aligned position ranges relative to the spliced CDS. Insertions appear as singletons. |
| Genome_length | Length of the ungapped genome sequence (including insertions) aligned via relaxed Smith-Waterman to reference. May be smaller than the original sequence file if divergent ends were hard-clipped. |
| Del_<AA/CDS/NT>_<start/end/len> | The start, end positions for amino acid, CDS, or genomic nucleotide deletions. Len is for total length. |
| Codon_shift | The number of extra nucleotide between the complete upstream codon and the insertion (0, 1 or 2). |
| In_frame | Specifies that the deletion contains no codon with partial deletions relative to `Ref ID`. |

## Special use of translated characters

Translation produces standard amino acid codes with the two non-standard exceptions listed below.  The translation engine also stops when it encounters a stop codon.

| Character | Interpretation |
| ------ | ------ |
| **.** | Missing alignment data (non-standard) |
| **-** | Gap in alignment (standard) |
| **~** | Partial codon (non-standard) |
| **X** | Ambiguous codon translation (standard) |

## CLI Usage

**Command-line** usage

```{bash}
Usage:

        ribosome                            help
                                            install (explicit)
                    [--module <MODULE>]     rebuild
                    [--module <MODULE>]     <†fasta|*tab> [<output1.seq> <output2.ins> <output3.del> [<output_file3.gen**] ]

        †if classified, fasta:  >ID|type_segment[_subtype]
        *if classified, tab:    ID<TAB>type_segment_[subtype]<TAB>sequence

        ** Also procduces output4.gen.ins and output5.gen.del if specified

        Valid modules: INFLUENZA, BETACORONAVIRUS

```

## Installation

For first time use, use the `install` command. Make sure you have at least *reporter* access to @vfn4's git repos:

- [convert](https://git.biotech.cdc.gov/vfn4/convert)
- [editMSA](https://git.biotech.cdc.gov/vfn4/editMSA)
- [sampling](https://git.biotech.cdc.gov/vfn4/sampling)
- [sswsort](https://git.biotech.cdc.gov/vfn4/sswsort)

Multiple references of the same type (and **same length**) can be included in the `spec/` subfolder.
References with new types may also be added, and must match reference types in `sswsort` for automatic classification purposes.
Automated classification is not needed if the reference type is included in the input (see below). Once references are updated, issue the `rebuild` command to deposit them for use in the `refs/` folder structure.

Influenza is the default module. Example usage might look like:

```{bash}
ribosome simple.fasta
ribosome flu.fasta out.seq out.ins
ribosome --module BETACORONAVIRUS cov.txt out.seq out.ins out.genome
```

Please note that the genome output is an additional optional argument given output files are explicitly specified (output files are normally sent to the calling working directory if not specified). An insertion file is implicitly created when genome output is specified. In this example it would be: `out.genome.ins`

## Updating the codon weight matrix

Given a curated alignment of some gene/genome and type, and, given the alignment has the same length as the references provided in the `spec/` and `refs/` for said gene/genome type, we can compute codon frequency statistics and store them to aid in codon re-adjustment / placement around *indels*. Given our repo is in `lib/`, one example could be:

```{bash}
lib/editMSA/storeStats.pl \
    <(himpala 'select reference_id,protein,cds_aln from sars_cov2.alignments limit 5000'| lib/convert/delim2fasta.pl) \
    -O BETACORONAVIRUS.sto -F 1,2
```

The more generic usage would be:

```{bash}
lib/editMSA/storeStats.pl \
    my-fasta-file-with-pipe-delimited-fields.fasta \
    -O BETACORONAVIRUS.sto -F 1,2
```

One can also display the codon weight matrix produced:

```{bash}
lib/editMSA/displayStats.pl  BETACORONAVIRUS.sto
```

## License

Generally this package is under GPLv3 but select files may be public domain. Datasets packaged with `dais-ribosome` are not for external re-distribution separate from this program. See the [LICENSE](LICENSE) file for more information.

## Component Citations

This package uses a modified version of SSW (`bin/licenses`):

> Zhao M, Lee WP, Garrison EP, Marth GT. SSW library: an SIMD Smith-Waterman C/C++ library for use in genomic applications. PLoS One. 2013;8(12):e82138. Published 2013 Dec 4. doi:10.1371/journal.pone.0082138

This package uses GNU Parallel :

> Tange O. GNU Parallel. Version 20200422. doi:10.5281/zenodo.1146014
