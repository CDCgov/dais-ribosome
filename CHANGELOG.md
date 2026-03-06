# DAIS-RIBOSOME Change Log

## v1.7.0 (TBD)

- **Changes**: The ranges for the `HA1` product for `A_HA_H1` and `A_HA_H3` are extended by 3 nucleotides and 18 nucleotides respectively

- **Fixes**:
  - The check for the required start codon in PB1-F2 now allows for AUG in addition to ATG
  - Fixes a bug causing some query coordinate ranges to be incorrectly shifted right for `M2`, `NEP`, and `PA-X`

## v1.6.2 (2026.03.19)

- **Fixes**: The check for the required start codon in PB1-F2 properly checks the beginning of the coding sequence instead of the query sequence

## v1.6.1 (2025.04.18)

- **Fixes**: The codon-weight matrix had the wrong number of sites for KF640637|L, this has been corrected.

## v1.6.0 (2025.04.18)

- **Feature**:
  - adds the A_NA_N4 reference, 'ONTARIO6118' to the INFLUENZA module.
  - adds the A_NA_N5 reference, 'RU1526' and 'ALASKA4733' (different lineage coordinate spaces) to the INFLUENZA module.
  - adds the A_NA_N6 reference, 'SICHUAN26221' to the INFLUENZA module.

## v1.5.6 (2025.03.28)

- **Changes**: ASTRAKHAN3212|A_NA_N8|OM403994 references now use lineage representatives LC339685 and MF973227, which enhances alignment quality.

## v1.5.5 (2024.10.30)

- **Changes**: Streamlines the dockerfile and changes to debian/bookworm.

## v1.5.4 (2024.08.22)

- **Change**: RSV ctypes are renamed, e.g., RSVA is RSV_AD; RSV reference IDs now match the NCBI accession.
- **Fixes**: updates the SC2 codon-weight matrix to induce corrected alignments for a deletion in S circa position 24.

## v1.5.3 (2024.07.25)

**Fixes**: SC2 alignments were sometimes missing biologically relevant indels. Smith-Waterman weights were relaxed to accommodate.

## v1.5.2 (2024.05.10)

- **Feature**:
  - adds docker support and adds usage
  - Adds configurable work directories using `IFX_WORK_DIR`.

- **Fixes**: fixes SC2 alignments for S protein insertions (MPLF) at position 16.

## v1.5.1 (2024.04.30)

**Fixes**: RSV annotations now correctly validated.

## v1.5 (2024.04.25)

**Feature**: adds the A_NA_N8 reference, 'ASTRAKHAN3212' to the INFLUENZA module.

## v1.4 (2024.04.03)

- **Feature**: adds built-in RSV module support. Many Thanks to C. Paden and P. Mandal!

## v1.3.3 (2023.06.26)

- **Change**: Right pads the CDS after premature stop codons.
- **Fixes**:
  - corrects coordinate map / range edge case involving insertions and incomplete CDS
  - "has insertion" boolean now reports `false` the only insertions are outside the terminated CDS
  - the CDS insertion table no longer provides insertions outside of the terminated CDS

## v1.3.2 (2023.06.23)

- Relaxes input file format detection and makes it invariant for all modules.

## v1.3.1 (2023.05.24)

- Restores right padding of CDS without adding back in the variable stop codon bugs from v1.2
- Premature stop codons will still be truncated.

## v1.3 (2023.05)

- Parallel SSW errors bubble up properly (including for fail-over from GRID jobs).
- Properly terminates CDS similar to AA sequences for premature stop codons.
- Skip empty exons and sequences during insertion processing.

## v1.2 (2022.08)

- Adds ORF7b and ORF9b to the BETACORONAVIRUS.spec (Thanks to C. Paden)

## v1.1 (2022.07)

- Cleaned up scripts and added individual licensing.

## v1.0 (2021.04)

- Initial tagged release supporting SARS-CoV-2 and Influenza.
