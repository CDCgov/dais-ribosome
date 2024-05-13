# DAIS-RIBOSOME Change Log

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
