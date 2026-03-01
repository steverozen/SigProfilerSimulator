# Indel Placement Algorithm in SigProfilerSimulator

## Where the mutation types and counts come from

The simulator does not invent mutation types or counts. They come from an **input catalog generated automatically from the user's real mutation data**:

1. The user provides real somatic mutations (VCF/MAF files) in the `input/` directory.
2. **SigProfilerMatrixGenerator** is called to classify those mutations and produce catalog matrix files — e.g., `project.ID83.all` for indels. These are TSV files where rows are mutation types (like `1:Del:T:5`) and columns are samples, with integer counts as values.
3. `mutation_preparation()` (in `mutational_simulator.py`) reads those catalog files into a dictionary: `samples[context][sample_name][mutation_type] = count`.
4. The simulator then generates exactly those counts of each type per sample, placing them at random matching positions.

This is what makes the tool a proper null hypothesis generator: it preserves "what kinds of mutations and how many" from the real data while randomizing "where they occur" on the genome.

## Chromosome allocation

Before placing individual mutations, `mut_tracker()` distributes indels across chromosomes **proportionally to chromosome size** (unlike SBS mutations which use nucleotide context distributions). This is because indels don't have a simple trinucleotide context — their context depends on repeat structure.

## Position selection: rejection sampling

The core algorithm is essentially **rejection sampling** — pick random positions until one matches the required sequence context:

1. **Pick a random position** uniformly on the chromosome using `fastrand.pcg32bounded(chromosome_length)`
2. **Extract the sequence** at that position (length = indel size)
3. **Skip** if the sequence contains N bases or overlaps an already-placed mutation
4. **Analyze the sequence context** at that position:
   - Count **tandem repeats** both forward and backward (up to 5 each direction)
   - Check for **microhomology** (partial sequence match at deletion boundaries)
   - Determine **subtype** from base composition: `T` (A/T bases), `C` (C/G bases), `R` (repeat), `M` (microhomology)
5. **Build a type key** like `"15T"` = 1bp indel, 5 repeats, T-type
6. **Check if this type matches** any remaining mutation in the catalog — if yes, place the mutation and write it out; if no, discard and loop back to step 1

## Three indel categories handled separately

The code splits indels into three groups, each with its own placement loop:

| Category | Condition | Example | Handling |
|---|---|---|---|
| **Standard** (`indel_types`) | All deletions + insertions with repeats > 0 | `1:Del:T:5` | Match position's repeat context to catalog |
| **Zero-repeat insertions** (`indel_types_O`) | Insertions with 0 repeats | `1:Ins:T:0` | Generate random inserted base(s); verify they differ from flanking sequence |
| **Microhomology** (`indel_types_M`) | Insertions at microhomology sites | `3:Ins:M:1` | Check for partial sequence match at boundaries |

## Repeat counting example

For a mutation like "1bp deletion of T in a homopolymer of length 5": the code picks a random position, checks if it's a T, then counts consecutive T's forward and backward. If `repeat_count_forward + repeat_count_reverse == 4` (the deleted T plus 4 copies = homopolymer of 5), the position matches and the deletion is placed there.

## Key implication

Because positions are sampled uniformly but **accepted only if the context matches**, the effective mutation density is higher in regions rich in the required context. For example, deletions requiring long homopolymers will cluster in homopolymer-rich regions — this naturally captures the biological tendency without explicitly modeling regional mutation rates.

## TSB-aware variant (ID-415)

The ID-415 classification extends this by partitioning positions by **transcriptional strand bias** (Transcribed, Untranscribed, Bidirectional, Non-genic). The random position is drawn only from genomic regions matching the required bias category, using precomputed bias region lengths and `bisect` for efficient lookup.

## Indel mutation type format

Indel mutations are represented as colon-separated strings:

```
[Length]:[Type]:[Subtype]:[Repeat_Count]
```

For example: `1:Del:T:5` = 1bp deletion, T-type base, 5 repeats

For ID-415 (with transcriptional strand bias), a bias prefix is added:

```
[Bias]:[Length]:[Type]:[Subtype]:[Repeat_Count]
```

Where Bias is one of: T (transcribed), U (untranscribed), B (bidirectional), N (non-genic), Q (mixed/unclassifiable).

## Key code locations

| Concept | File | Lines (approx) |
|---|---|---|
| Indel type parsing | `mutational_simulator.py` | 1248-1262 |
| Three-category split | `mutational_simulator.py` | 1266-1295 |
| Random position selection | `mutational_simulator.py` | 1308-1323 |
| Repeat counting (forward/reverse) | `mutational_simulator.py` | 1352-1406 |
| Subtype determination | `mutational_simulator.py` | 1408-1457 |
| Deletion matching & output | `mutational_simulator.py` | 1465-1625 |
| Insertion matching & output | `mutational_simulator.py` | 1627-1826 |
| Microhomology handling | `mutational_simulator.py` | 1828-2343 |
| Zero-repeat insertion loop | `mutational_simulator.py` | 2620-2882 |
| ID-415 TSB simulation | `mutational_simulator.py` | 2885-4600 |
| Chromosome allocation | `mutational_simulator.py` | 843-915 |
