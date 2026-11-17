## Tax-credit Step (reference database benchmarking)

The other three steps process your data. The tax-credit step answers a question that comes
*before* that: **which reference database, classify method and confidence threshold should I
use for this marker gene?**

It benchmarks any number of reference databases against any number of classify methods and
parameter sets, scores the results, and produces summary tables and plots. It does not need
outputs from the qaqc, repseqs or taxonomy steps.

### Why benchmark

Taxonomic assignment is the step where results most depend on choices you cannot check by
eye. A classifier will happily assign a species name that the reference database could never
have gotten right, and a confidence threshold that is well calibrated for 16S may be far too
permissive for 12S. Benchmarking tells you, for your marker and your candidate databases, how
often assignments are correct, how often they are too shallow, and how often they are simply
wrong.

### Evaluation methods

| Method | What it does | Answers |
|---|---|---|
| `cross-validated` | Taxonomy-aware K-fold splits of a database; classify held-out sequences. | How well does this database classify sequences like the ones it contains? |
| `cross-validated-trad` | Random splits of the query list only; the reference keeps every sequence. | A near-best-case ceiling on a random subset (see the caveat below). |
| `novel-taxa` | Hold out whole taxa, so the query's own taxon is absent from the reference. | What happens to organisms your database has never seen? |
| `self-validated` | Classify the full database against itself. | Best case ceiling; catches internal inconsistencies. |
| `mock-community` | Classify real sequencing data from communities of known composition. | How does it do on real reads, including PCR and abundance effects? |

The first four are simulated from the reference database itself. `mock-community` needs real
data you supply: a feature table, ASV sequences, and the expected composition and/or the known
taxonomy of each ASV.

> **`cross-validated-trad` does not hold sequences out of the reference.** Each fold's
> `ref_seqs.fasta` and `ref_taxa.tsv` are symlinks to the full database, so every query is
> classified against a reference that still contains it, exact self-match included. Only the
> query *list* changes between folds, and a single classifier is fitted per database and reused
> across all of them. Read its scores as a ceiling on a random subset rather than as
> cross-validated performance, and use `cross-validated` when you want a genuine held-out split.

By default every fold of `cross-validated`, `cross-validated-trad` and `novel-taxa` scores all of
its queries: the database is divided between the `iterations` folds, about
`n_sequences / iterations` sequences each. On a large database such as SILVA that is the
dominant cost of the step. `max_query_size_per_fold` caps it: any fold holding more queries than
the cap scores a reproducible random subset of that size instead. It is a per-fold count, not a
total. Only the queries shrink: references are unchanged, and queries left out are never added
back to them. Novel-taxa folds are built from the full cross-validated folds before capping. The
uncapped set stays in each capped fold as `query_full.fasta` / `query_full_taxa.tsv`, and
`self-validated` is never capped. It replaces `trad_cv_query_size`, which a config can no longer
set. Changing it on a finished run needs the datasets phase to run again (delete
`.datasets.done`). See [Configuration](../configuration.md) for the details.

### Requirements

The step calls the sibling [tour2-tax-credit](https://github.com/ksil-NOAA/tour2-tax-credit) package. Clone
its `tour2-tax-credit` branch; there is nothing to install, because Tourmaline imports it directly
from `tax_credit_package_dir`:

```bash
cd ..    # alongside your tourmaline clone
git clone -b tour2-tax-credit https://github.com/ksil-NOAA/tour2-tax-credit.git tax-credit
```

Cloning into a directory named `tax-credit` keeps the shipped configs working as-is, since
`tax_credit_package_dir` defaults to `../tax-credit`. Any other name works too — just point
`tax_credit_package_dir` at it (relative paths are resolved from the Tourmaline directory). Keep
the clone up to date with Tourmaline; the two change together.

The classify methods
use the same environments as the taxonomy step, so `bt2-blca` and `revamp` need theirs — see
[Install and Setup](../install.md).

### Running

```bash
conda activate snakemake-tour2
./tourmaline.sh -s tax-credit -c config_04_tax_credit.yaml -n 8
```

Benchmark runs are much larger than a normal analysis — a full matrix of databases × methods ×
parameter sets × folds is hundreds of assignment jobs. Start with a reduced matrix.

> **Keep `classify_threads` well below `--cores`.** Assignment jobs run in parallel, and a
> large `classify_threads` lets a single job claim every core, serializing the run.

Not every stage uses the threads it reserves, so `tax_credit_assign_fold` sizes each job
individually rather than giving every job `classify_threads`:

| job | threads |
|---|---|
| `classify-sklearn`, `classify-consensus-blast`, `classify-consensus-vsearch`, bowtie2 index + alignment, the REVAMP BLAST | `classify_threads` |
| `makeblastdb`, the BLCA stage, the confidence reformat, REVAMP's post-BLAST assignment | 1 |
| `fit-classifier-naive-bayes`, as a shared fit job (confidence sweeps, `cross-validated-trad`) | 1 |

`fit-classifier-naive-bayes` takes no thread option. A fold with a single confidence fits and
classifies in one job, which gets `classify_threads` for the classify step.

### Memory: you must opt in, or nothing is throttled

Because single-threaded jobs ask for one core, **threads no longer limit how many run at
once — memory does.** These jobs are not comparable in size: the BLCA stage loads the
whole fold reference into memory before running muscle, and a naive-bayes fit over a full
database with a wide `fit_params` ngram-range is larger still. Sixteen of either on a
16-core box will exhaust the machine.

The rule therefore reserves memory per job type, from the `mem_mb_*` config keys. But
**Snakemake only enforces a resource when the run declares its ceiling**, so pass `--mem`:

```bash
./tourmaline.sh -s tax-credit -c config_04_tax_credit.yaml -n 16 --mem 120000
```

Set it to roughly the RAM you are willing to give the run (`free -m`). Without `--mem` the
reservations are recorded and ignored, and the run can OOM. The wrapper always passes
`--retries 2`, which helps because each attempt multiplies a job's reservation, so a job
killed for memory gets more on the next try.

Driving Snakemake directly, the equivalent is:

```bash
snakemake --use-conda -s tax_credit_step.Snakefile run_tax_credit \
  --configfile config_04_tax_credit.yaml --cores 16 --latency-wait 15 \
  --resources mem_mb=120000 --retries 2
```

The shipped defaults are starting points, not measurements — memory scales with the
reference. Every job's peak memory is recorded automatically (see
[Benchmarks](#benchmarks-where-the-time-goes)): after a run, read `peak_rss_mb` (or
`peak_pss_mb` on Linux, the better figure for classify-sklearn's worker processes) for each
job kind in `summaries/benchmark_summary.tsv`, and set the matching `mem_mb_*` key a little
above it.

> **`--resources mem_mb` must be at least as large as the biggest single job.** A run
> whose ceiling is below some job's reservation fails at execution with that job, and a
> dry run will not catch it: Snakemake evaluates neither threads nor resources under `-n`.

When a run dies during the fitting or BLCA phase, the signature is an OOM one level down
rather than a clean error — `muscle: error while loading shared libraries: libc.so.6:
cannot map zero-fill pages`, `*** OUT OF MEMORY ***` after only tens of MB, then
`ValueError: No records found in handle` where BLCA read muscle's empty output.

### consensus-blast builds a BLAST database per fold

`classify-consensus-blast` only honours `--p-num-threads` against a pre-indexed database.
Given `--i-reference-reads` it falls back to `blastn -subject`, which is single-threaded and
warns that the thread count is ignored. Tourmaline therefore runs `makeblastdb` for each fold
and classifies with `--i-blastdb`.

None of the swept consensus parameters change the database, so a sweep builds it once per fold
in a fit job that every parameter combination then reuses — the same sharing bt2-blca uses for
its bowtie2 index. A fold with a single parameter combination builds it in-job instead. The
databases land in `data/results/<evaluation-method>/<fold>/<reference>/consensus-blast/consensus_blast/blastdb.qza`
and are roughly the size of the reference, so a large sweep over many folds needs the disk for
one copy per fold.

> **Scores shift slightly versus runs made before this change.** BLAST computes E-values from
> the effective size of the search space, and `-db` and `-subject` mode size it differently, so
> a borderline hit can fall on either side of the `evalue` cutoff. Re-run consensus-blast
> rather than comparing new numbers against old ones.

### Smoke tests

Small subsampled MiFish databases in `00-data/tax-credit-test/` exercise every code path
quickly. They are **test fixtures only** — never report benchmarking results from them.

```bash
./tourmaline.sh -s tax-credit -c config_04_tax_credit_test.yaml -n 8       # reduced matrix
./tourmaline.sh -s tax-credit -c config_04_tax_credit_test_full.yaml -n 8  # full matrix, slower
./tourmaline.sh -s tax-credit -c config_04_tax_credit_test_mock.yaml -n 8  # mock-community only
```

See `00-data/tax-credit-test/README.md` for what a complete run looks like and how the
fixtures are regenerated.

> On these tiny databases, novel-taxa scores collapse toward zero. That is an expected
> artifact of subsampling, not a bug: removing a taxon from a small reference often removes
> its whole lineage.

### REVAMP is mock-community only

`revamp` classifies against NCBI `nt`, which is far too large to simulate cross-validated,
novel-taxa or self-validated datasets from. It therefore runs **only** for `mock-community`.
Listing it alongside other evaluation methods is fine — it is skipped for them with a notice.
Listing it with no `mock-community` in `evaluation_methods` stops the run with an error.

### Gotcha: staging does not re-run when only the config changes

The evaluation plan is written by the staging phase, which takes the config as parameters
rather than as an input file. If you **add a mock dataset or reference database** to an
existing run name, delete the `.datasets.done` marker in the run output directory, or:

```bash
snakemake --use-conda -s tax_credit_step.Snakefile run_tax_credit \
  --configfile config_04_tax_credit.yaml --cores 8 \
  --forcerun tax_credit_prepare_datasets
```

Tourmaline checks for this and stops with an explanatory message rather than silently
skipping the new database.

### Running one phase directly

For debugging, individual phases can be run without Snakemake:

```bash
python scripts/run_tax_credit.py --config config_04_tax_credit.yaml --phase evaluate
```

Phases are `datasets`, `manifest`, `evaluate`, `plot`, `post-assign` and
`benchmark-summary`.

### Benchmarks: where the time goes

Every rule in the step records a Snakemake benchmark file for each job: wall time, CPU time,
peak memory and mean CPU load, measured over the job's whole process tree (`qiime`,
bowtie2, muscle, classify-sklearn's workers). At the end of the run,
`tax_credit_benchmark_summary` joins them to the manifest and writes two tables to
`summaries/`:

| File | Contents |
|---|---|
| `benchmark_summary.tsv` | Totals per evaluation method × classify method × job kind, plus one row per phase (dataset generation, manifest, mock scoring, evaluation, plotting). Sorted by `total_wall_s`, so the first row is the dominant cost. |
| `job_benchmarks.tsv` | One row per assignment job: `job_id`, evaluation method, database, dataset, classify method, `job_kind`, and the measurements. |

`job_kind` says what a job actually ran: `fit` (classifier, bowtie2 index or BLAST database
only), `classify` (against a shared fit), `fit+classify` (both, in one job),
`blca-align` / `blca-blca` / `blca-reformat` (the three bt2-blca stages; `blca-blca` is
muscle plus bootstraps), or `blca-all-stages` (a single-combination bt2-blca fold).

The columns:

- `total_wall_s`, `median_wall_s`, `max_wall_s`: wall time in seconds.
- `total_cpu_s`: CPU seconds.
- `peak_rss_mb` / `peak_pss_mb`: the largest single job's peak memory in MB. PSS is Linux only.
- `mean_load_pct`: average CPU use, where 100 = one full core. Compare it with the threads a
  job reserved: a single-threaded BLCA job reads about 100, and a threaded job reading
  far below `classify_threads × 100` is not using its threads.
- `pct_of_total_wall`: each row's share of all measured wall time.

The raw per-job files stay in `benchmarks/` (`assign/<job_id>.tsv`,
`mock_score/<job_id>.tsv`, and one file per phase).

Caveats:

- **On macOS only wall time is recorded.** Snakemake reads processes through psutil's
  `memory_full_info()`, which macOS refuses even for your own processes, so it writes 0 for
  CPU time, memory and I/O. The summary detects this and leaves those columns blank rather
  than reporting 0 MB. Linux, including HPC nodes and Docker, reports every column, so size
  `mem_mb_*` from a Linux run.
- **Failed and killed jobs write nothing**, and a retried job records only the attempt that
  succeeded. A job killed for memory twice before succeeding looks like one run.
- **Memory is sampled**, every 0.5 s for a job's first ~15 s and every 30 s after that, so a
  short spike in a long job can be missed.
- **Wall time includes waiting** for CPU or disk when jobs compete. Compare `total_cpu_s`
  with `total_wall_s`.
- **Runs from before benchmarking was added** are not rerun to fill the files in. Snakemake
  only regenerates a missing benchmark file for its target rule, so their jobs appear with
  blank measurements (counted by `n_measured`), and the summary rule alone runs.

The cost is negligible: one small file per job, and a sampling thread that reads OS counters.

### Outputs

```
[run_name]-tax-credit/
├── data/           # staged reference databases, simulated datasets, assignment results
├── summaries/      # metric tables, per-job scores, best-run selections, benchmark summaries
├── benchmarks/     # raw per-job Snakemake benchmark files (time, CPU, peak memory)
└── plots/          # one folder per evaluation method
```

Each assignment result directory under `data/` also holds `classification_accuracy_log.tsv`
(one line per query) and `classification_metrics.tsv` (one row per rank).

### How assignments are scored

Every query, or every mock ASV, is a *match*, an *underclassification* (correct but too
shallow), an *overclassification* (deeper than the reference knows) or a *misclassification*.
Precision, Recall and F-measure are computed from those tallies, with one change from the
original tax-credit: **overclassification is a false positive but no longer a false negative**,
so it lowers precision and F-measure but not recall. The original scoring is still reported,
as `Recall_original` and `F-measure_original`, for comparison; `match_ratio` equals
`Recall_original`. Scores that have nothing to score are NaN rather than 0.

Results scored before this change are not comparable with new ones. Re-score a finished run
with `--forcerun tax_credit_evaluate`; no classification is repeated. See
[Configuration](../configuration.md), *Scoring* and *Mock community evaluation*, for every
column and the full definitions.

### Full configuration reference

Every key — reference databases, evaluation methods, simulation settings, parameter sweeps,
mock-community inputs, metrics and plotting — is documented in
[Configuration](../configuration.md), section *4. Tax-credit configuration*.
