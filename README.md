
The output for the `dais-ribosome` consists of two tab-delimited files. One with `.ins` for insertions and `.seq` for sequence related data.

The insertion file example:<br />

| ID | C_type | Ref_ID | Protein | Site | Codon | Residue |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 11209 | B_HA | PHUKET3073 | HA | 161 | aaa | K | 
| 154957 | B_HA | PHUKET3073 | HA | 163 | krc | X | 
| 223550 | B_HA | PHUKET3073 | HA | 161 | caa | Q | 

The sequence file header:

| ID | C_type | Ref_ID | Protein | VH | Insertion | AA_seq | AA_aln | cds_id | NT_seq | NT_aln | Query_nt_coordinates | CDS_nt_coordinates |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 223550 | B_HA | BRISBANE60 | HA-signal | e81d2d895c70e91bb3ef917fe49fdab7 | false | MKAIIVLLMVVTSNA | MKAIIVLLMVVTSNA | 2aa6443b92ca45b301faa4d46e5fbd3b010e3ab7 | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCA | 20..64 | 1..45 |
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
| NT_seq / NT_aln | The nucleotide cds sequence (less deletions + insertions) and the cds alignment (bases relative to reference only). |
| Query_nt_coordinates | Set of aligned position ranges representing the aligned coordinates relative to the original submitted query sequence. Insertions appear as singletons. |
| CDS_nt_coordinates | Set of aligned position ranges relative to the spliced CDS. Insertions appear as singletons. |