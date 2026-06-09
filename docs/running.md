## Running

The main entrypoint is the `tourmaline.sh` script. Run from the `tourmaline` directory.

### Usage

```bash
conda activate snakemake-tour2
./tourmaline.sh --step [qaqc,repseqs,taxonomy,tax-credit] --configfile [config1,config2,...] --cores N
```

Notes:

- The number of `--step` entries must match the number of `--configfile` entries.
- Provide steps and configs in the same order.
- Each Snakemake rule uses `--use-conda` to pull the `qiime2-amplicon-2024.10` env as needed.

### Examples

Run a single step (taxonomy):

```bash
./tourmaline.sh -s taxonomy -c config_03_taxonomy.yaml -n 6
```

Run all steps with one command:

```bash
./tourmaline.sh -s qaqc,repseqs,taxonomy -c config_01_qaqc.yaml,config_02_repseqs.yaml,config_03_taxonomy.yaml -n 6
```

Run tax-credit only (reference database benchmarking):

```bash
./tourmaline.sh -s tax-credit -c config_04_tax_credit.yaml -n 6
```

### Tips

- Use `--printshellcmds` and `--dryrun` when calling Snakemake directly to see planned commands and DAG.
- Choose rarefaction depth after inspecting table summaries in Repseqs step.
- Optionally filter by taxonomy/length/abundance/prevalence before final visualizations.

### Output structure

```
v2-results/
├── [run_name]-qaqc/
├── [run_name]-repseqs/
├── [run_name]-taxonomy/
└── [run_name]-tax-credit/
```

See [Configuration](configuration.md) for parameters and [Steps](steps/qaqc.md) for per-step details.


