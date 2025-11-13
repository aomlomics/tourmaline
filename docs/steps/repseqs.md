## Repseqs (ASV Generation)

Generates ASVs using DADA2 (paired/single) or Deblur (single), with optional filtering and optional diversity visualizations; produces feature tables, representative sequences, and supporting summaries.

### Key configuration

```yaml
asv_method: dada2pe  # dada2pe | dada2se | deblur

# DADA2 examples
dada2_trunc_len_f: 0
dada2pe_trunc_len_r: 0
dada2_trim_left_f: 0
dada2pe_trim_left_r: 0

# Deblur examples
deblur_trim_length: 150
deblur_trim_left: 0

# Filtering (optional)
to_filter: false
repseq_min_length: 0
repseq_max_length: 100000
repseq_min_abundance: 0.0
repseq_min_prevalence: 0.0
```

Input sources:

- From QA/QC: reuse the same `run_name`/`output_dir` or set `qaqc_run_name` to point at another QA/QC run.
- External demultiplexed reads: provide `fastq_qza_file` (QIIME 2 artifact).
- Metadata: set `sample_metadata_file` to a TSV/CSV (auto-generated if omitted).

Optional diversity outputs (alpha rarefaction + core metrics) require

```yaml
plot_diversity: true
alpha_max_depth: 500
core_sampling_depth: 500
```

### Outputs

- `[run_name]-table.qza`, `[run_name]-repseqs.qza`
- Stats/visualizations under `stats/` including denoising summaries, feature-table summaries, sequence-length reports, and diversity plots (if enabled)

Proceed to [Taxonomy](taxonomy.md).


