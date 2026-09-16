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

The `NA`-bearing lineages include internal `NA` ranks (e.g.
`Eukaryota;Chordata;Actinopteri;NA;Lutjanidae;Lutjanus;Lutjanus griseus`), which
makes these databases a useful check of rank-position handling — see
[NA ranks](#na-ranks) below.

## Running

```bash
conda activate snakemake-tour2
export R_HOME="$CONDA_PREFIX/../qiime2-amplicon-2024.10/lib/R"   # see note below

./tourmaline.sh -s tax-credit -c config_04_tax_credit_test.yaml -n 8       # reduced matrix, ~30 min
./tourmaline.sh -s tax-credit -c config_04_tax_credit_test_full.yaml -n 8  # full matrix, slower
```

Outputs land in `../v2-results/tax_credit_test-tax-credit/` and
`../v2-results/tax_credit_test_full-tax-credit/`.

## What a complete reduced run looks like

With `config_04_tax_credit_test.yaml` as of 2026-09-15 (all three evaluation
methods, 2 CV folds, novel-taxa levels 5 and 6, two values per parameter sweep):

| Check | Expected |
|---|---|
| Snakemake log | `145 of 145 steps (100%) done` |
| Datasets | 14: 2 self-validated, 4 cross-validated, 8 novel-taxa |
| `assignment_manifest.tsv` / `assignment-done/` | 140 jobs = 14 datasets × 10 (2 naive-bayes, 2 consensus-vsearch, 4 bt2-blca, plus a naive-bayes fit and a bowtie2 index per dataset) |
| `summaries/` | `evaluate_classification_summary_{CV,novel,self_validated}.csv`, every database × method × parameter set present, no NaN F-measures |
| `plots/` | 25 files each for cross-validated, novel-taxa, self-validated |
| bt2-blca `taxonomy.tsv` | one row per query in `staging/query.fasta` |

Most of this can be checked from the log, the manifest and directory listings,
without opening the per-read result files. `log_analysis/novel-taxa`
`method_parameter_sensitivity-species.csv` and `-genus.csv` are identical by
design: novel-taxa expected lineages stop at genus or family.

**Runs from before 2026-09-15 have wrong bt2-blca metrics and wrong internal-NA
handling in `log_analysis/`** (both fixed, see below). Rerun the evaluate, plot
and log-analysis phases — and reformat the bt2-blca `raw-taxonomy.tsv` files —
before using older outputs.

## Gotchas this test surfaced

1. **`novel_taxa_levels` and `novel_recall_min_level` must agree.** Dataset
   generation uses `novel_taxa_levels`; the assignment manifest independently
   enumerates `range(6, novel_recall_min_level, -1)` and never checks the fold
   directories exist. A mismatch fails with `ref_seqs.qza does not exist`.
   `novel_taxa_levels: [6]` needs `novel_recall_min_level: 5`; `[5,6]` needs `4`.

2. **Keep `classify_threads` well below `--cores`.** At `classify_threads: 4`
   with `-n 6`, two jobs would need 8 threads, so snakemake runs assignment jobs
   one at a time — the same work took 1022 s serialized versus 462 s at
   `classify_threads: 2` with `-n 8`.

3. **`self-validated` + `bt2-blca` used to crash (fixed).** Self-validation makes
   the query set identical to the reference set, so a query could reach scoring
   with no usable hits and `scripts/blca_from_bowtie.py` called `max()` on an
   empty score dict. The script now writes `Unclassified` for such queries, and
   the reduced config runs `self-validated` with bt2-blca.

4. **bt2-blca confidence truncation left empty ranks (fixed).**
   `scripts/reformat_summary_for_r.py` used to keep every rank whose own
   confidence passed the cutoff and pad the rest with empty strings, producing
   `…;Hexagrammidae;;` or internal gaps like `…;Actinopteri;;Centropomidae;…`.
   tax-credit scores those as misclassifications, which is why raising the BLCA
   confidence appeared to *increase* misclassification. It now stops at the
   first rank below the cutoff and writes only the kept ranks
   (`…;Hexagrammidae`). A correct bt2-blca `taxonomy.tsv` has no empty ranks.

## NA ranks

Reference lineages can carry `NA` at any rank. tax-credit treats them as follows
(`tax_credit.taxa_manipulator.normalize_taxon`):

- **Trailing** `NA` or empty ranks mean "not assigned below here" and are
  dropped: `…;Lutjanus;NA` is a genus-level lineage.
- **Internal** `NA` ranks are kept so every rank stays at its position:
  `Actinopteri;NA;Lutjanidae` is still class / order / family. Never strip them —
  doing so shifts genus into the family slot and moves errors to the wrong level.

Before 2026-09-15 tax-credit removed every `;NA`, so `classification_accuracy_log.tsv`
was correct but `log_analysis/` tables printed `…;Actinopteri;Lutjanidae;…`,
reported the species as `expected_genus`, and per-level ratios counted genus
errors at family level.

Unit tests for this live in the tax-credit package (`tax_credit/tests/`, see the
tourmaline `CLAUDE.md` for the command).

## R_HOME note

On a machine with a system R framework of a different architecture, `qiime`
fails to start (`incompatible architecture … libR.dylib`) because rpy2 resolves
`R_HOME` to `/Library/Frameworks/R.framework` instead of the conda env's R.
Export `R_HOME=<qiime env>/lib/R` before running.
