The DAIS **ribosome** compartmentalizes the original translation engine developed for our protein analytics database. Currently there is only support for INFLUENZA, but some seeds for extensibility have been planted throughout the code. I provide a brief outline of the algorithm:
1.  If necessary, classify the nucleotide gene segment into its influenza type, segment, and subtype (compound or `C_type`) as well as make sure CDS is on the *forward strand*†.
2.  Align said segments (via [SSW](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082138)) to the corresponding reference sequence(s) and pick the best alignment.
3.  There can be more than one reference reading frame per C_type, so complete step 2 for each `reference_id`. Data with no reference is held-aside and added in later.
4.  Fix alignment ends that have been chopped (due to local alignment disagreement).
5.  Create protein product CDS using an internal specification.
5.  Correct product alignment so that indels occur within frame only, then tabularize and create `cds_id`.
6.  Amend insertion tables to use protein coordinates and translate.
7.  Translate CDS to amino acids, calculating the `variant_hash` as well.
8.  Create coordinate mapping between CDS and AA
9.  Combine AA, CDS, and coordinate tabular data; output with insertion data from step 6.

***

**Input** can be one of the four formats.

1.  Unannotated Fasta (ID only)
>\>223550<br />ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA
2.  Annotated Fasta (ID and compound type)
>\>223550|B_HA<br />ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA
3.  Unannotated tab-delimited (ID only)
>223550 *\<tab\>* ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA
4.  Annotated tab-delimited (ID and compound type)
>223550 *\<tab\>* B_HA *\<tab\>* ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA

†Annotated fasta / tab-delimited data **must** be on the forward or plus strand. Unannotated data will be classified and reverse complemented to the forward strand in the CDS as needed.
***

**Output** for the `dais-ribosome` consists of two tab-delimited files. One with `.ins` for insertions and `.seq` for sequence related data. An insertion file output example:<br />

| ID | C_type | Ref_ID | Protein | Site | Codon | Residue |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 11209 | B_HA | PHUKET3073 | HA | 161 | aaa | K | 
| 154957 | B_HA | PHUKET3073 | HA | 163 | krc | X | 
| 223550 | B_HA | PHUKET3073 | HA | 161 | caa | Q | 

A sequence file output example:

| ID | C_type | Ref_ID | Protein | VH | Insertion | AA_seq | AA_aln | CDS_id | CDS_seq | CDS_aln | Query_nt_coordinates | CDS_nt_coordinates |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 223550 | B_HA | BRISBANE60 | HA-signal | e81d2e81d2d895c70e91bb3ef917fe49fdab7d89549fdab7 | false | MKAIIVLLMVVTSNA | MKAIIVLLMVVTSNA | 2aa6443b92ca45b301faa4d46e5fbd3b010e3ab7 | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | 20..64 | 1..45 |
| 223550 | B_HA | PHUKET3073 | HA-signal | e81d2d895c70e91bb3ef917fe49fdab7 | false | MKAIIVLLMVVTSNA | MKAIIVLLMVVTSNA | 2aa6443b92ca45b301faa4d46e5fbd3b010e3ab7 | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | 20..64 | 1..45 |
| 11209 | B_HA | BRISBANE60 | HA-signal | c7ee7ff234abf5c0591e0fe1af26ca87 | false | MKAIIILLMVVTSNA | MKAIIILLMVVTSNA | c49a73ab7280362c8c710abbf648708c41f97712 | ATGAAGGCAATAATTATACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTATACTACTCATGGTAGTAACATCCAATGCA | 1..45 | 1..45 |
| 11209 | B_HA | PHUKET3073 | HA-signal | c7ee7ff234abf5c0591e0fe1af26ca87 | false | MKAIIILLMVVTSNA | MKAIIILLMVVTSNA | c49a73ab7280362c8c710abbf648708c41f97712 | ATGAAGGCAATAATTATACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTATACTACTCATGGTAGTAACATCCAATGCA | 1..45 | 1..45 |

The field explanations:

| field | description |
| ------ | ------ |
| ID | Any unique identifier, but likely the `flu_sequence_id`, `epi_segment_id`, raw `nt_id`, or NCBI accession. |
| C_type | The compound type conisting of the influenza type, segment, and subtype if applicable. This is the same as IRMA. Chimeric types start with an asterisk. | 
| Ref_ID | As with DAIS, the reference reading frame used for alignment. |
| Protein | The protein product or peptide derived from the gene segment. |
| VH | The `variant_hash` as used in DAIS (md5 hex of `AA_seq`). |
| Insertion | Boolean indicating whether or not there is an insertion relative to reference. |
| Site | The upstream amino acid position relative the insertion. |
| Codon | The codon(s) inserted. |
| Residue | the residue(s) inserted. Partial codons are denoted as `?`. |
| AA_seq / AA_aln | The amino acid sequence (less deletions + insertions) and the amino acid alignment (residues relative to reference only). |
| CDS_seq / CDS_aln | The nucleotide cds sequence (less deletions + insertions) and the cds alignment (bases relative to reference only). |
| Query_nt_coordinates | Set of aligned position ranges representing the aligned coordinates relative to the original submitted query sequence. Insertions appear as singletons. |
| CDS_nt_coordinates | Set of aligned position ranges relative to the spliced CDS. Insertions appear as singletons. |