## Running

The main entrypoint is the `tourmaline.sh` script. Run it from the Tourmaline repository root —
rules invoke `scripts/...` by relative path. Outputs can go anywhere, via `output_dir`.

### Usage

```bash
conda activate snakemake-tour2
./tourmaline.sh --step [qaqc,repseqs,taxonomy] --configfile [config1,config2,...] --cores N
```

Short flags: `-s` / `-c` / `-n`.

Notes:

- The number of `--step` entries must match the number of `--configfile` entries.
- Provide steps and configs in the same order.
- Each Snakemake call uses `--use-conda` to pull the `qiime2-amplicon-2024.10` env as needed,
  plus `--latency-wait 15` for shared filesystems.

### Examples

Run a single step (taxonomy):

```bash
./tourmaline.sh -s taxonomy -c config_03_taxonomy.yaml -n 6
```

Run all steps with one command:

```bash
./tourmaline.sh -s qaqc,repseqs,taxonomy -c config_01_qaqc.yaml,config_02_repseqs.yaml,config_03_taxonomy.yaml -n 6
```

### Running Snakemake directly

Call Snakemake yourself for dry runs, single rules, or `--printshellcmds`. Each step has its own
Snakefile and target rule:

| Step | Snakefile | Target rule |
|---|---|---|
| qaqc | `qaqc_step.Snakefile` | `qaqc_all` |
| repseqs | `repseqs_step.Snakefile` | `run_denoise` |
| taxonomy | `taxonomy_step.Snakefile` | `run_taxonomy` (the default target) |

```bash
snakemake --use-conda -s qaqc_step.Snakefile     qaqc_all     --configfile config_01_qaqc.yaml --cores 6 --dryrun
snakemake --use-conda -s repseqs_step.Snakefile  run_denoise  --configfile config_02_repseqs.yaml --cores 6
snakemake --use-conda -s taxonomy_step.Snakefile run_taxonomy --configfile config_03_taxonomy.yaml --cores 6
```

You can also request a single output file as the target, which is the quickest way to run one
optional rule — for example the per-base quality drop-off report, which is not part of
`qaqc_all`:

```bash
snakemake --use-conda -s qaqc_step.Snakefile \
  ../v2-results/my_run-qaqc/stats/my_run-seq_qual_dropoff.txt \
  --configfile config_01_qaqc.yaml --cores 4
```

### Reusing upstream runs

Steps chain through the filesystem, not through Snakemake, so a later step reads whatever `.qza`
an earlier one wrote. Input resolution follows a fixed precedence:

1. An explicit upstream run name (`qaqc_run_name`, `repseqs_run_name`) — read from that run's directory.
2. An explicit external artifact (`fastq_qza_file`; `repseqs_qza_file` + `table_qza_file`).
3. Otherwise, this step's own `run_name` directory.

This is what lets you re-run one step over several parameter sets while sharing upstream output —
for example, three taxonomy runs with different classifiers, all pointing at one `repseqs_run_name`.

### Parameter sweeps

To compare many parameter sets, expand a base config over a parameter space, then run the
generated configs in parallel:

```bash
python scripts/generate_configs.py <base_config> <parameter_space_config>

scripts/run_parallel_tourmaline.sh \
  --config-dir parameter_sweep_configs \
  --config-prefix config-01-qaqc \
  --step qaqc \
  --parallel-jobs 4 \
  --cores-per-job 6
```

Example parameter space files are in `00-data/` (`parameter_space.yaml`,
`parameter_space_repseqs.yaml`, `parameter_space_taxonomy.yaml`).

The parameter space file controls the naming of what is generated:

| Key | Default | Effect |
|---|---|---|
| `output_dir` | `parameter_sweep_configs` | Where the generated configs are written |
| `config_prefix` | `config-01-qaqc` | Generated files are `{config_prefix}_000.yaml`, `_001`, … |
| `run_name_prefix` | `test_data` | Each config gets `run_name: {run_name_prefix}_000`, … |

Pass the same `config_prefix` to `run_parallel_tourmaline.sh` as `--config-prefix`, and point
`--config-dir` at the same `output_dir`. Because each config gets its own `run_name`, the results
land in separate output directories and can be compared directly.

Parallel execution requires GNU `parallel`, which is included in the `snakemake-tour2`
environment. Total cores used is roughly `--parallel-jobs × --cores-per-job`, so size it to the
machine.

To compare *classifiers or databases* rather than denoising parameters, the run-name chaining
above is usually simpler: point several taxonomy configs at one `repseqs_run_name`.

### HPC / SLURM

`scripts/sbatch_tourmaline2_step*.sh` are SLURM array wrappers for the same steps. Adapt the
account, partition, time and memory directives to your cluster before submitting.

### Tips

- Use `--dryrun` before any long run; use `--printshellcmds` to see the exact commands.
- Choose rarefaction depth after inspecting the table summaries from the Repseqs step.
- Give each attempt a distinct `run_name` so results sit side by side instead of overwriting.
- Rules have no `log:` directives — output goes to the terminal. Redirect it if you want a record.

### Output structure

```
output_dir/
├── [run_name]-qaqc/
├── [run_name]-repseqs/
└── [run_name]-taxonomy/
```

Each step copies its config file into its own output directory as `{run_name}-{step}_config.yaml`,
so a run's provenance lives next to its results. This is also what
[`format_analysisMetadata.py`](metadata.md) reads.

See [Configuration](configuration.md) for parameters and [Steps](steps/qaqc.md) for per-step details.
