## Flu Module Notes

### References & Flu Subtypes

The flu module currently supports influenza A and B. 

For influenza B, most of the flu segments use the reference PHUKET3073 which represents B/Yamagata, although the HA and NS segments also use BRISBANE60 which is a B/Victoria reference. 

For influenza A, all segments other than the hemagglutinin and neuraminidase use references from HK4801. The following references handle the other hemagglutinin subtypes:

- H1: CALI07 (H1N1)
- H2: ANNARBOR60 (H2N2)
- H3: HK4801 (H3N2)
- H5: VT1203
- H7: ANHUI01 (H7N9)
- H9: BGD0994

The other neuraminidase subtypes use the references:

- N1: CALI07 (H1N1)
- N2: HK4801 (H3N2)
- N4: ONTARIO6118
- N5: RU1526 and ALASKA4733
- N6: SICHUAN26221
- N7: NL219
- N8: ASTRAKHAN3212

### Protein Products

All eight flu segments are included as compound types by the flu module. They generate the following products:

- Segment 1 (`PB2`): The whole protein is represented as a product `PB2`. `PB2-S1` is not currently included.
- Segment 2 (`PB1`): The whole protein is represented as a product `PB1`. Through leaky scanning, `PB1-F2` can be produced, which is a virulence factor. `PB1-N40` is also made via leaky scanning but is not currently included.
- Segment 3 (`PA`): The whole protein is represented as a product `PA`. In influenza A specifically, `PA-X` is produced via a +1 ribosomal frameshift, which plays a role in host shutoff. `PA-N182` and `PA-N155` (with truncated N-terminus) are not currently included.
- Segment 4 (`HA`): `HA-signal` is the signal peptide, followed by `HA` (the full HA protein before cleavage, after the signal peptide has been removed) and `HA1` (the globular head, post-cleavage). `HA1` is represented as a distinct product because it is used downstream in analyses such as protein modeling, and its epitopes are of interest. `HA2` (the stalk) is not currently included.
- Segment 5 (`NP`): The whole protein is represented as a product `NP`.
- Segment 6 (`NA`): The whole protein is represented as a product `NA`.
- Segment 7 (`MP`): The `M1` matrix protein is included. For influenza A, the alternatively-spliced membrane protein `M2` is also included. For influenza B, the membrane protein `BM2` is produced via a stop-start pentanucleotide. The rare alternative splicing `M42` for influenza A is not included.
- Segment 8 (`NS`): The nonstructural protein `NS1` and alternatively-spliced nuclear export protein `NEP` are both included.

### Stop Codons

In previous releases of DAIS-ribosome, there was inconsistency with whether the `HA` products should include the stop codon or not. The flu module now includes stop codons for consistency, and to enable stop extension.
