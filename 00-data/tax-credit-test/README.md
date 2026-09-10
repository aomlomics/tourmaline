# tax-credit test databases

Small, subsampled MiFish reference databases for smoke-testing the tax-credit
step. **These are test fixtures, not production databases** — do not use them to
produce real benchmarking results.

| File | Sequences | Source database |
|---|---|---|
| `mifishGom-test-{seqs,taxa}.{qza,fasta,tsv}` | 400 | `gomex-mifish-trimAJM-250815-derepU-seqs.qza` + `gomex-mifish-trimAJM-260106-derepU-taxa.qza` (2,093 seqs) |
| `addJ-test-{seqs,taxa}.{qza,fasta,tsv}` | 500 | `addJ-rc_12S_efc-fish-addSI-addMF-seqs-derepU.qza` + `…-taxa-derepU.qza` (16,954 seqs) |

Both `.qza` and plain FASTA/TSV forms are provided; `run_tax_credit.py` accepts
either (a `.qza` is exported to text on the fly, a text path is used as-is).

## Regenerating

Deterministic — the same seed reproduces these files byte for byte:

```bash
export R_HOME="$CONDA_PREFIX/lib/R"   # only if a host R framework shadows the env's

conda run -n qiime2-amplicon-2024.10 python scripts/subsample_reference_db.py \
    --seqs <full-db>-seqs.qza --taxa <full-db>-taxa.qza \
    --n-target 400 --n-multi-species 50 \
    --out-prefix 00-data/tax-credit-test/mifishGom-test --write-qza

conda run -n qiime2-amplicon-2024.10 python scripts/subsample_reference_db.py \
    --seqs <full-db>-seqs.qza --taxa <full-db>-taxa.qza \
    --n-target 500 --n-multi-species 60 \
    --out-prefix 00-data/tax-credit-test/addJ-test --write-qza
```

Seed `20260910` (the script default). See the module docstring in
[`scripts/subsample_reference_db.py`](../../scripts/subsample_reference_db.py)
for why the sample is stratified rather than random: `get_strata()` needs
species carrying at least `iterations` sequences to form meaningful
cross-validation strata, and novel-taxa evaluation needs sibling taxa to remain
in the reference set when a genus/family/order is held out.

## Structure retained

| | mifishGom-test | addJ-test |
|---|---|---|
| sequences | 400 | 500 |
| species (with ≥3 seqs) | 268 (50) | 347 (60) |
| genera | 251 | 322 |
| families (with ≥2 genera) | 169 (54) | 215 (73) |
| orders (with ≥2 families) | 58 (35) | 66 (48) |
| `NA`-bearing lineages | 20 | 56 |

Sequence strings, FASTA headers (including the
`_representative_of_N_identical_accessions` suffix convention) and 7-rank
taxonomy strings are copied verbatim from the source databases, in source file
order, so these are drop-in replacements for the full-size files.

## Running

```bash
conda activate snakemake-tour2
export R_HOME="$CONDA_PREFIX/../qiime2-amplicon-2024.10/lib/R"   # see note below

./tourmaline.sh -s tax-credit -c config_04_tax_credit_test.yaml -n 8       # reduced, ~8 min
./tourmaline.sh -s tax-credit -c config_04_tax_credit_test_full.yaml -n 8  # full matrix, slower
```

Outputs land in `../v2-results/tax_credit_test-tax-credit/` and
`../v2-results/tax_credit_test_full-tax-credit/`.

Verified run of the reduced config: 29/29 snakemake steps, 24/24 assignment
jobs, 462 s wall clock, producing CV and novel-taxa summaries plus 50 plots for
both databases across all three classify methods.

## Three gotchas this test surfaced

1. **`novel_taxa_levels` and `novel_recall_min_level` must agree.** Dataset
   generation uses `novel_taxa_levels`; the assignment manifest independently
   enumerates `range(6, novel_recall_min_level, -1)` and never checks the fold
   directories exist. A mismatch fails with `ref_seqs.qza does not exist`.
   For `novel_taxa_levels: [6]` the matching value is `novel_recall_min_level: 5`.

2. **Keep `classify_threads` well below `--cores`.** At `classify_threads: 4`
   with `-n 6`, two jobs would need 8 threads, so snakemake runs assignment jobs
   one at a time — the same work took 1022 s serialized versus 462 s at
   `classify_threads: 2` with `-n 8`.

3. **`self-validated` + `bt2-blca` crashes.** Self-validation makes the query set
   identical to the reference set, so a query whose only bowtie2 hit is itself
   yields an empty score dict and `scripts/blca_from_bowtie.py` calls `max()` on
   it unguarded (`ValueError: max() arg is an empty sequence`). Reachable on any
   database; near-certain on a sparse one. The reduced config omits
   `self-validated` for this reason.

## R_HOME note

On a machine with a system R framework of a different architecture, `qiime`
fails to start (`incompatible architecture … libR.dylib`) because rpy2 resolves
`R_HOME` to `/Library/Frameworks/R.framework` instead of the conda env's R.
Export `R_HOME=<qiime env>/lib/R` before running.
