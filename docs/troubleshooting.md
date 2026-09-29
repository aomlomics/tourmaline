## Troubleshooting and FAQ

### Configuration errors

**`KeyError: 'some_key'` before anything runs.**
The Snakefiles read most config values directly, so a *missing* key fails at parse time. Copy a
fresh example config and edit it, rather than writing one from scratch or deleting lines. An
optional key should be present but empty:

```yaml
repseqs_qza_file:        # correct — present, empty, means "not set"
```

**A parameter change had no effect.**
Check that you edited the config you actually passed to `--configfile`. Each run also copies its
config into its output directory as `{run_name}-{step}_config.yaml` — compare that against what
you expected.

**Steps and config files are mismatched.**
`tourmaline.sh` requires the same number of `--step` and `--configfile` entries, in the same
order. The script stops with a usage message if they differ.

### Environment errors

**A rule fails saying the conda environment cannot be found.**
Rules request environments **by name**, so they must already exist on the machine. The names
matter exactly: `qiime2-amplicon-2024.10`, and `bt2-blca`, `revamp` or `krona` for those
features. Check with `conda env list`.

**`snakemake: command not found`.**
Activate the Snakemake environment first: `conda activate snakemake-tour2`.

**A rule cannot find `scripts/...`.**
Rules invoke scripts by relative path, so you must run from the Tourmaline repository root.
Output can still go anywhere via `output_dir`.

**R fails to load `dplyr` or `Biostrings` (taxonomy step).**
A host R installation on `PATH` can shadow the conda environment's R. Export `R_HOME` to point
at the QIIME 2 environment's R:

```bash
export R_HOME="$(conda info --base)/envs/qiime2-amplicon-2024.10/lib/R"
```

### Input and data errors

**Manifest paths are wrong.** Manifest files need *absolute* paths, and the headers must match
the expected format exactly. See [QA/QC step](steps/qaqc.md#manifest-formats).

**Samples are missing after import.** Without a manifest, files must follow the expected naming:
`{sample}_R1.fastq.gz` / `{sample}_R2.fastq.gz`, or the `_R1_001.fastq.gz` variant. The sample
name is the file name prefix.

**Deblur and sample names.** Avoid underscores in sample names; use alphanumerics, dashes and
periods.

**A QIIME 2 artifact will not load.** Artifacts written by one QIIME 2 version are often **not**
readable by earlier versions. Import and run with 2024.10.

### Results look wrong

**Filtering removed every sequence.** If `to_filter: True`, `repseq_max_length` must be a real
upper bound — the filter keeps `length <= repseq_max_length`, so the template default of `0`
removes everything. See [Configuration](configuration.md).

**Almost all reads are lost at the DADA2 step.** Usually the truncation lengths are too
aggressive, leaving paired reads unable to overlap. Check `stats/dada2_stats.tsv` to see which
stage lost the reads, and revisit the Step 1 quality plots. Forward + reverse truncation lengths
must exceed the amplicon length by enough to overlap (roughly 12 bp minimum, plus a margin).

**Everything is `Unassigned`.** Check that the reference database matches your marker gene and
that `taxa_ranks` matches the database's rank structure. For naive-bayes, a classifier trained
on a different primer region will perform poorly. Lower `skl_confidence` only after ruling
these out.

**REVAMP `mostEnvOUT` results look like `allIN`.** Taxid filtering needs BLAST's `taxdb` files
on the `BLASTDB` path. Without them BLAST prints a warning and **continues with the exclusion
list unapplied**. Tourmaline checks for this before starting; see
[Taxonomy step](steps/taxonomy.md#revamp).

### Tax-credit

**Jobs run one at a time.** Keep `classify_threads` well below `--cores`, or a single
assignment job claims every core and the rest serialize.

**A newly added database or mock dataset produces no jobs.** Staging does not re-run when only
the config changes. Delete the `.datasets.done` marker in the run output directory, or pass
`--forcerun tax_credit_prepare_datasets`. See [Tax-credit step](steps/tax_credit.md).

**Novel-taxa scores are near zero on the test databases.** Expected. The subsampled fixtures in
`00-data/tax-credit-test/` are too small for meaningful novel-taxa results — they are for
smoke-testing only.

### Debugging technique

- **Dry run first.** `--dryrun` shows what Snakemake plans to do without doing it.
- **See the actual commands.** `--printshellcmds` prints each shell command. Rules have no
  `log:` directives, so output goes to the terminal.
- **Re-run one rule.** Call Snakemake directly with the step's Snakefile and a single target:

  ```bash
  snakemake --use-conda -s repseqs_step.Snakefile run_denoise \
    --configfile config_02_repseqs.yaml --cores 6 --printshellcmds --dryrun
  ```

- **Force regeneration.** Delete the specific output file, or use `--forcerun <rule>`.
- **Compare runs.** Change `run_name` rather than overwriting, so both results remain.

### Performance

- Set `--cores` to match the machine. Denoising and classification parallelize well.
- Test truncation and trimming parameters on a subset before committing to a full dataset.
- For many parameter combinations, use the sweep tooling in [Running](running.md) rather than
  editing configs by hand.

### Where to ask questions

- Tourmaline GitHub issues: <https://github.com/aomlomics/tourmaline/issues>
- [QIIME 2 Forum](https://forum.qiime2.org) for questions about the underlying QIIME 2 commands
