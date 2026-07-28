<!-- omit from toc -->
# DAIS-ribosome

DAIS-ribosome or **ribosome** annotates CDS and protein products for supported virus genomes into database-oriented output.

**Table of contents:**

- [Input data](#input-data)
- [Output data](#output-data)
- [Field descriptions](#field-descriptions)
- [Special use of translated characters](#special-use-of-translated-characters)
- [CLI Usage](#cli-usage)
  - [Bare Metal](#bare-metal)
  - [Container](#container)
    - [Verifying images on `ghcr.io`](#verifying-images-on-ghcrio)
  - [Grid](#grid)
- [Installation](#installation)
- [Methodology](#methodology)
- [Notices](#notices)
  - [Contact Info](#contact-info)
  - [Public Domain Standard Notice](#public-domain-standard-notice)
  - [License Standard Notice](#license-standard-notice)
  - [Privacy Standard Notice](#privacy-standard-notice)
  - [Contributing Standard Notice](#contributing-standard-notice)
  - [Records Management Standard Notice](#records-management-standard-notice)
- [Additional Standard Notices](#additional-standard-notices)

## Input data

**Input** can be one of the four formats.

1. Unannotated FASTA (ID only)

    ```fasta
    >223550
    ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA
    ```

2. Annotated FASTA (ID and compound type)

    ```fasta
    >223550|B_HA
    ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA
    ```

3. Unannotated tab-delimited (ID only)

    ```text
    223550<TAB>ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA
    ```

4. Annotated tab-delimited (ID and compound type)

    ```text
    223550<TAB>B_HA<TAB>ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA
    ```

†Annotated FASTA / tab-delimited data **must** be on the forward or plus strand. Unannotated data will be classified and reverse-complemented to the forward strand as needed.

## Output data

**Output** for `dais-ribosome` consists of three product tab-delimited files: `.seq.txt` for sequence-related data, `.ins.txt` for insertions, and `.del.txt` for deletions. Genome `.gen_seq.txt`, `.gen_ins.txt`, and `.gen_del.txt` files are also written unless skipped. An insertion file output example:

| query_id | ctype | reference_id | product_name | upstream_aa_pos | inserted_nt | inserted_aa | upstream_nt_pos | codon_shift |
| -------- | ----- | ------------ | ------------ | --------------- | ----------- | ----------- | --------------- | ----------- |
| 11209    | B_HA  | PHUKET3073   | HA           | 161             | AAA         | K           | 483             | 0           |
| 154957   | B_HA  | PHUKET3073   | HA           | 163             | KRC         | X           | 489             | 0           |
| 223550   | B_HA  | PHUKET3073   | HA           | 161             | CAA         | Q           | 483             | 0           |

A deletion file example:

| query_id       | ctype      | reference_id | product_name | variant_hash                     | del_aa_start | del_aa_end | del_aa_len | in_frame | cds_id                                   | del_cds_start | del_cds_end | del_cds_len |
| -------------- | ---------- | ------------ | ------------ | -------------------------------- | ------------ | ---------- | ---------- | -------- | ---------------------------------------- | ------------- | ----------- | ----------- |
| EPI_ISL_410721 | SARS-CoV-2 | WUHAN19      | orf1ab       | 5ba70e95c9a3251bc6155f62295dd3e8 | 994          | 1002       | 9          | true     | 29cd767e2d144c31179395fd606d1489ce731746 | 2980          | 3006        | 27          |
| EPI_ISL_410721 | SARS-CoV-2 | WUHAN19      | orf1ab       | 5ba70e95c9a3251bc6155f62295dd3e8 | 1012         | 1012       | 1          | true     | 29cd767e2d144c31179395fd606d1489ce731746 | 3034          | 3036        | 3           |
| EPI_ISL_410721 | SARS-CoV-2 | WUHAN19      | S            | 450c068c437e7536d27fdb883d95d4f4 | 72           | 72         | 1          | true     | 36a75a0d34960c048abaf82ee46a1b713eee534e | 214           | 216         | 3           |
| EPI_ISL_410721 | SARS-CoV-2 | WUHAN19      | S            | 450c068c437e7536d27fdb883d95d4f4 | 146          | 146        | 1          | true     | 36a75a0d34960c048abaf82ee46a1b713eee534e | 436           | 438         | 3           |
| EPI_ISL_410721 | SARS-CoV-2 | WUHAN19      | S            | 450c068c437e7536d27fdb883d95d4f4 | 254          | 256        | 3          | true     | 36a75a0d34960c048abaf82ee46a1b713eee534e | 760           | 768         | 9           |
| EPI_ISL_410721 | SARS-CoV-2 | WUHAN19      | S            | 450c068c437e7536d27fdb883d95d4f4 | 680          | 683        | 4          | true     | 36a75a0d34960c048abaf82ee46a1b713eee534e | 2038          | 2049        | 12          |

A sequence file output example:

| query_id | ctype | reference_id | product_name | variant_hash                     | aa_seq          | aa_aln          | cds_id                                   | has_insertion | has_shift_indel | cds_seq                                       | cds_aln                                       | query_coordinates | cds_coordinates |
| -------- | ----- | ------------ | ------------ | -------------------------------- | --------------- | --------------- | ---------------------------------------- | ------------- | --------------- | --------------------------------------------- | --------------------------------------------- | ----------------- | --------------- |
| 223550   | B_HA  | BRISBANE60   | HA-signal    | e81d2d895c70e91bb3ef917fe49fdab7 | MKAIIVLLMVVTSNA | MKAIIVLLMVVTSNA | 2aa6443b92ca45b301faa4d46e5fbd3b010e3ab7 | false         | false           | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | 20..64            | 1..45           |
| 223550   | B_HA  | PHUKET3073   | HA-signal    | e81d2d895c70e91bb3ef917fe49fdab7 | MKAIIVLLMVVTSNA | MKAIIVLLMVVTSNA | 2aa6443b92ca45b301faa4d46e5fbd3b010e3ab7 | false         | false           | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | 20..64            | 1..45           |
| 11209    | B_HA  | BRISBANE60   | HA-signal    | c7ee7ff234abf5c0591e0fe1af26ca87 | MKAIIILLMVVTSNA | MKAIIILLMVVTSNA | c49a73ab7280362c8c710abbf648708c41f97712 | false         | false           | ATGAAGGCAATAATTATACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTATACTACTCATGGTAGTAACATCCAATGCA | 1..45             | 1..45           |
| 11209    | B_HA  | PHUKET3073   | HA-signal    | c7ee7ff234abf5c0591e0fe1af26ca87 | MKAIIILLMVVTSNA | MKAIIILLMVVTSNA | c49a73ab7280362c8c710abbf648708c41f97712 | false         | false           | ATGAAGGCAATAATTATACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTATACTACTCATGGTAGTAACATCCAATGCA | 1..45             | 1..45           |

Genome output file example:

| query_id      | ctype      | reference_id | genome_id                                | genome_length | has_insertion | genome_seq                                                                                                       | genome_aln                                                                                           |
| ------------- | ---------- | ------------ | ---------------------------------------- | ------------- | ------------- | ---------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| EPI_ISL_FAKE1 | SARS-CoV-2 | WUHAN19      | 8e193a72b22a666947b21cb785af6780c2c6996b | 108           | true          | TTTAAGGTTTATACCTTCCCAGGTAACAAACCAACC**TGGGTTTGG**AACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTC | .TTTAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTC |

Genome insertion output file example:

| query_id      | ctype      | reference_id | upstream_nt_pos | inserted_nt |
| ------------- | ---------- | ------------ | --------------- | ----------- |
| EPI_ISL_FAKE1 | SARS-CoV-2 | WUHAN19      | 37              | TGGGTTTGG   |

Genome deletion output file example:

| query_id      | ctype      | reference_id | del_start | del_end | del_len |
| ------------- | ---------- | ------------ | --------- | ------- | ------- |
| EPI_ISL_FAKE1 | SARS-CoV-2 | WUHAN19      | 3246      | 3272    | 27      |

## Field descriptions

| field                             | description                                                                                                                                                                                        |
| --------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| query_id                          | Any unique identifier for the input sequences, but likely the `flu_sequence_id`, `epi_segment_id`, raw `nt_id`, or NCBI accession.                                                                 |
| ctype                             | The compound type, which represents the classification level used by the module. This is the type/segment/subtype for `flu` and the taxon for `cov`.                                               |
| reference_id                      | The ID for the coordinate reference group which was aligned against.                                                                                                                               |
| protein                           | The protein product or peptide derived from the gene segment.                                                                                                                                      |
| variant_hash (aa_id)              | The amino acid sequence ID using the MD5 hash of `aa_seq`.                                                                                                                                         |
| has_insertion                     | Boolean indicating whether there is an insertion relative to reference in the original CDS.                                                                                                        |
| has_shift_indel                   | Boolean indicating whether an insertion or deletion in the CDS is a non-triplet and could have induced a frameshift.                                                                               |
| upstream_aa_pos / upstream_nt_pos | The upstream amino acid / nucleotide position for the insertion relative to the reference coordinates.                                                                                             |
| inserted_aa / inserted_nt         | The amino acids / nucleotides inserted.                                                                                                                                                            |
| cds_id / genome_id                | The nucleotide sequence ID using the sha1 hex of `cds_seq` and `genome_seq`.                                                                                                                       |
| aa_seq / aa_aln                   | The amino acid sequence (less deletions + insertions) and the amino acid alignment (residues relative to reference only).                                                                          |
| cds_seq / cds_aln                 | The nucleotide CDS sequence (less deletions + insertions) and the CDS alignment (bases relative to reference only).                                                                                |
| genome_seq / genome_aln           | The nucleotide genome sequence (less deletions + insertions) and the genome alignment (bases relative to reference only).                                                                          |
| query_coordinates                 | Set of position ranges in the original submitted query sequence used to form `cds_seq`, including inserted query segments.                                                                         |
| cds_coordinates                   | Set of position ranges relative to the spliced CDS. Insertions appear as singleton insertion positions.                                                                                            |
| genome_length                     | Length of the ungapped genome sequence (including insertions) aligned via relaxed Smith-Waterman to reference. May be smaller than the original sequence file if divergent ends were hard-clipped. |
| del_<aa/cds/nt>_<start/end/len>   | The start, end positions for amino acid, CDS, or genomic nucleotide deletions. Len is for total length.                                                                                            |
| codon_shift                       | The number of extra nucleotides between the complete upstream codon and the insertion (0, 1, or 2).                                                                                                |
| in_frame                          | Specifies that the deletion contains no codon with partial deletions relative to `reference_id`.                                                                                                   |

## Special use of translated characters

Translation produces [standard amino acid codes] with the two non-standard
exceptions listed below.  The translation engine also stops when it encounters a
stop codon.

| Character | Interpretation                                                                                                                                                                  |
| --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **.**     | On left: padding to ensure aligned sequence outputs are in the proper coordinates. On right (optional): padding to ensure the output spans the full reference coordinate space. |
| **-**     | Gap in alignment                                                                                                                                                                |
| **~**     | Partial codon                                                                                                                                                                   |
| **X**     | Ambiguous codon translation                                                                                                                                                     |
| **\***    | Stop codon                                                                                                                                                                      |

_Note:_ `.` conveys custom semantics; `~` is non-IUPAC.

## CLI Usage

### Bare Metal

Execute `./ribosome --help` for the latest options. For example:

```bash
Usage: ribosome [OPTIONS] <DATA_FILE> [SEQUENCE_OUTPUT] [INSERTION_OUTPUT] [DELETION_OUTPUT] [GENOME_SEQ] [GENOME_INS] [GENOME_DEL]

Arguments:
  <DATA_FILE>
          Data file to annotate in TSV or FASTA format.†

  [SEQUENCE_OUTPUT]
          CDS and AA output, including coordinate mapping information, as a filename or path

  [INSERTION_OUTPUT]
          Insertion output filename or path

  [DELETION_OUTPUT]
          Deletion output filename or path

  [GENOME_SEQ] [GENOME_INS] [GENOME_DEL]
          Genome sequence, insertion, and deletion output paths. Passing a single genome output prefix still works but is deprecated

Options:
      --output-prefix <OUTPUT_PREFIX>
          The prefix to use for naming the output files (or an existing folder in which to place them)

  -m, --module <MODULE>
          Name of the alignment module

          [default: flu]

  -T, --threads <THREADS>
          Run in simultaneous multi-threaded mode

  -G, --is-grid-task
          Automatically detect the array task id from SGE environment variables and write partition files for downstream collation.

          Output files are required and will be suffixed with a partition id.

  -S, --submit-grid-job <SUBMIT_GRID_JOB>
          Submit and block on a grid engine (SGE) array job of the specified size

      --verbose
          Prints warning messages to stderr

      --assume-default-ctype <ASSUME_DEFAULT_CTYPE>
          A default ctype to use if any input records are not annotated. If not specified, an SSWSort module will be used to classify the query if a module exists, otherwise an error is produced
```

### Container

Images are available both on ghcr.io and Docker hub.

```bash
# Use /tmp in the container
docker run --rm -v $(pwd):/data -t cdcgov/dais-ribosome:latest ribosome flu.fasta t1.seq.txt t1.ins.txt t1.del.txt t1.gen

# Alter the scratch directory to use our mount
docker run --rm -v $(pwd):/data -e IFX_WORK_DIR=/data -t cdcgov/dais-ribosome:latest ribosome flu.fasta
```

#### Verifying images on `ghcr.io`

While we publish to both `ghcr.io` and Docker Hub, please use the former for
cryptographic verification of the image signature and SLSA build provenance
attestation. Verification with Sigstore [cosign](https://github.com/sigstore/cosign):

```bash
# Replace with the version of interest
TAG=test

cosign verify --new-bundle-format ghcr.io/cdcgov/dais-ribosome:$TAG \
  --certificate-oidc-issuer https://token.actions.githubusercontent.com \
  --certificate-identity-regexp '^https://github[.]com/CDCgov/dais-ribosome/[.]github/workflows/release[.]yml@refs/tags/.+$'

cosign verify-attestation --new-bundle-format --type slsaprovenance1 \
  ghcr.io/cdcgov/dais-ribosome:$TAG \
  --certificate-oidc-issuer https://token.actions.githubusercontent.com \
  --certificate-identity-regexp '^https://github[.]com/CDCgov/dais-ribosome/[.]github/workflows/release[.]yml@refs/tags/.+$'
```

### Grid

***DAIS-ribosome*** can submit jobs on SGE for you and block on the
results. Output partitions are concatenated together. One must ensure a network
file system is accessible from each node, otherwise the user can write a more
sophisticated wrapper using `--is-grid-task` for the executors.

```bash
# 100 tasks or partitions are created for the job array but results will be
# concatenated and removed at the end.
./ribosome flu.fasta out.seq.txt out.ins.txt out.del.txt --submit-grid-job 100
```

## Installation

***DAIS-ribosome*** [is packaged in all releases since v2.0.0][releases]. A release can be
used with the following steps:

1. Download the `tar.gz` file for your platform from the [release
   page][releases]
2. Unzip it with `tar -xzfv path/to/tar.gz`
3. The binary can be run within that folder, such as `./ribosome --help`. If the
   binary is moved to a new location, ensure that `ribosome_res` and
   `sswsort_res` are also moved.

Alternatively, you can compile and install DAIS-ribosome yourself with the
following steps:

1. Install [rustup]
2. Add the `nightly` toolchain with `rustup toolchain install nightly`
3. Clone DAIS-ribosome using `git clone https://github.com/CDCgov/dais-ribosome`
4. Enter the directory for DAIS-ribosome with `cd dais-ribosome`
5. Run `./install.sh` to compile DAIS-ribosome and fetch SSWSort resources

For RHEL 8 compatible Linux distributions, you can either rebuild from the
included
[`Dockerfile`](https://github.com/CDCgov/dais-ribosome/blob/main/Dockerfile) or
pull our latest pre-built image:

```bash
# Since v2.0.0
docker pull ghcr.io/cdcgov/dais-ribosome:latest

# Includes legacy v1 versions
docker pull docker.io/cdcgov/dais-ribosome:latest
```

## Methodology

I provide a brief outline of the algorithm:

1. If necessary, classify the nucleotide gene segment / genome into a module-recognized label. For influenza this is a compound type (or `ctype`) made up of flu type, segment, and subtype; it could also be different species or variants as the module allows. If a classifier assigns a reverse-strand match, the query is reverse-complemented; otherwise the input is expected to already be on the *forward strand*†.
2. For each coordinate `reference_id` available for the label, align the query to the corresponding reference sequence(s) and pick the best alignment. Currently [Zoe]'s striped Smith-Waterman implementation is used.
3. Apply module rules around alignment, such as trimming to a plausible start, repairing locally clipped ends, or extending to an in-frame stop codon.
4. Convert the genomic alignment states into query/reference coordinate ranges for matches, deletions, and insertions.
5. Intersect those ranges with the configured exon/CDS specification for each protein product, producing product-level query/CDS coordinate ranges.
6. Correct eligible product indels by shifting coordinate ranges to in-frame codon boundaries, then merge remaining adjacent deletions and skip products missing a required start codon.
7. Materialize the product ranges into CDS and amino acid sequences, alignments, coordinate mappings, insertions, and deletions. Product sequences stop at the first in-frame stop codon outside insertions; insertion rows filter singletons, doubletons, and entirely ambiguous inserts.
8. Produce product CDS alignment, insertion, and deletion tables, and optionally produce genomic alignment, insertion, and deletion tables.

## Notices

### Contact Info

For direct correspondence on the project, feel free to contact: [Samuel S. Shepard](mailto:sshepard@cdc.gov), Influenza Division, National Center for Immunization and Respiratory Diseases, Centers for Disease Control and Prevention or reach out to other [contributors](CONTRIBUTORS.md).

### Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).  All contributions to this repository will be released under the CC0 dedication.  By submitting a pull request you are agreeing to comply with this waiver of copyright interest.

### License Standard Notice

The repository utilizes code licensed under the terms of the Apache Software License and therefore is licensed under ASL v2 or later. This source code in this repository is free: you can redistribute it and/or modify it under the terms of the Apache Software License version 2, or (at your option) any later version. This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache Software License for more details. You should have received a copy of the Apache Software License along with this program. If not, see: <http://www.apache.org/licenses/LICENSE-2.0.html>. The source code forked from other open source projects will inherit its license.

### Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and information. All material and community participation is covered by the [Disclaimer](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md). For more information about CDC's privacy policy, please visit <http://www.cdc.gov/other/privacy.html>.

### Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo) and submitting a pull request. (If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git).) By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the [Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or later.

All comments, messages, pull requests, and other submissions received through CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Records Management Standard Notice

This repository is not a source of government records, but is a copy to increase collaboration and collaborative potential. All government records will be published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).

[standard amino acid codes]: https://www.bioinformatics.org/sms/iupac.html
[Zoe]: https://github.com/CDCgov/zoe
[releases]: https://github.com/CDCgov/dais-ribosome/releases
[rustup]: https://forge.rust-lang.org/infra/other-installation-methods.html