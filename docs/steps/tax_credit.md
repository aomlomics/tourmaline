## Tax-credit Step (reference database benchmarking)

> **In development.** The tax-credit step currently lives on the `feature/tax-credit-module`
> branch and is not yet part of the main `V2` branch.

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
| `cross-validated-trad` | Traditional random K-fold splits. | The same, without taxonomy-aware stratification. |
| `novel-taxa` | Hold out whole taxa, so the query's own taxon is absent from the reference. | What happens to organisms your database has never seen? |
| `self-validated` | Classify the full database against itself. | Best case ceiling; catches internal inconsistencies. |
| `mock-community` | Classify real sequencing data from communities of known composition. | How does it do on real reads, including PCR and abundance effects? |

The first four are simulated from the reference database itself. `mock-community` needs real
data you supply: a feature table, ASV sequences, and the expected composition and/or the known
taxonomy of each ASV.

### Requirements

The step calls the sibling [tax-credit](https://github.com/aomlomics/tax-credit) package, which
must be installed into the QIIME 2 environment:

```bash
conda activate qiime2-amplicon-2024.10
pip install -e ../tax-credit
```

Point `tax_credit_package_dir` at that clone (default `../tax-credit`). The classify methods
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

Phases are `datasets`, `manifest`, `evaluate`, `plot` and `post-assign`.

### Outputs

```
[run_name]-tax-credit/
├── data/           # staged reference databases, simulated datasets, assignment results
├── summaries/      # metric tables, per-job scores, best-run selections
└── plots/          # one folder per evaluation method
```

### Full configuration reference

Every key — reference databases, evaluation methods, simulation settings, parameter sweeps,
mock-community inputs, metrics and plotting — is documented in
[Configuration](../configuration.md), section *4. Tax-credit configuration*.
