# Summary: Bergstrom et al. (2020) — Generating realistic null hypothesis of cancer mutational landscapes using SigProfilerSimulator

**Citation:** Bergstrom EN, Barnes M, Martincorena I, Alexandrov LB. *BMC Bioinformatics* 21:438 (2020). https://doi.org/10.1186/s12859-020-03772-3

## Overview

SigProfilerSimulator is a computationally efficient tool for generating realistic simulated somatic mutational landscapes of cancer genomes. It constructs null hypothesis models by randomly shuffling observed mutations while preserving their mutational signatures at a user-selected resolution. This enables downstream statistical testing and evaluation of bioinformatics tools.

## How it works

- Takes real somatic mutations (VCF/MAF) and produces simulated mutations while maintaining mutational burden and mutational pattern at a preselected resolution.
- Uses a Monte Carlo approach with a PCG random number generator.
- Mutations are projected onto chromosomes proportionally based on the observed rate of each mutational channel in a reference genome.
- Works with SigProfilerMatrixGenerator to classify mutations before simulation.

## Mutation classification resolutions

| Classification | Context | Channels |
|---|---|---|
| SBS-6 | Mutation type only | 6 |
| SBS-24 | + transcriptional strand bias (TSB) | 24 |
| SBS-96 | + 1 bp flanking context | 96 |
| SBS-384 | + 1 bp + TSB | 384 |
| SBS-1536 | + 2 bp flanking context | 1,536 |
| SBS-6144 | + 2 bp + TSB | 6,144 |
| SBS-24576 | + 3 bp flanking context | 24,576 |

Also supports indel classifications: ID-83 and ID-415.

## Key simulation options

1. **Gender-aware** — proper handling of sex chromosomes
2. **Transcriptional strand bias** — accounts for transcription-coupled NER
3. **Updating mode** — mutations treated as dependent sequential events (each mutation updates trinucleotide counts in the reference)
4. **Chromosome-level preservation** — maintain per-chromosome mutational burden/pattern
5. **Exome simulations** — mutations only in protein-coding regions
6. **Poisson noise** — add noise to mutation counts per channel
7. **Probability mask** — increase/decrease mutation opportunity in specific genomic regions
8. **Germline simulation** — simulate variants in matched-normal samples

## Three application results

### 1. Doublet base substitutions (DBSs) are single genomic events
- Simulated 2,144 PCAWG genomes 1,000 times each at SBS-96 resolution.
- Compared expected adjacent SBS pairs (by chance) to observed DBSs.
- Found 10- to 1,000-fold enrichment of real DBSs over simulated.
- Conclusion: the vast majority of DBSs are not due to two adjacent SBSs occurring by chance; they are likely single genomic events.

### 2. Extended sequence context (+/-2 bp) is needed for mutational signatures
- Simulated PCAWG at SBS-6, SBS-96, and SBS-1536; compared to real data at higher resolutions.
- SBS-6 simulations fail to capture +/-2 bp patterns (91% of samples below 0.85 cosine similarity).
- SBS-96 simulations partially capture +/-2 bp (44% below 0.85).
- SBS-1536 captures +/-3 bp patterns well (only 6.5% below 0.85).
- Conclusion: SBS-1536 (+/-2 bp) is the sweet spot; extending to +/-3 bp (SBS-24576) adds little.

### 3. False-positive rates of driver gene detection tools
- Simulated 1,024 TCGA breast cancer WES samples 100 times at SBS-384 and SBS-6144.
- Tested MutSigCV1.41, MutSigCV2, and dNdScv on simulated data (which should have no drivers).
- MutSigCV1.41: 1.3-1.6 false-positive driver genes per simulation (q < 0.10).
- MutSigCV2: 0.2-0.3 false-positive driver genes per simulation.
- dNdScv: 0.02-0.03 false-positive driver genes per simulation (best performance).
- Lowering threshold to q < 0.01 eliminates all false positives for dNdScv and MutSigCV2.

## Performance

- Simulated ~37 million somatic mutations from 2,144 PCAWG genomes in ~90 seconds.
- Benchmarked on dual Intel Xeon Gold 6132 (2.60 GHz), 192 GB RAM.

## Availability

- Python: https://github.com/AlexandrovLab/SigProfilerSimulator
- R wrapper: https://github.com/AlexandrovLab/SigProfilerSimulatorR
- Docs: https://osf.io/usxjz/wiki/home/
- License: BSD 2-Clause
