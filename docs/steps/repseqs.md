## Repseqs (ASV Generation)

Generates ASVs using DADA2 (paired/single) or Deblur (single), with optional filtering; produces feature table and representative sequences.

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

Input source:

- From QA/QC: use the same `run_name`/`output_dir` or set `sample_run_name`
- External: set `fastq_qza_file` to a demultiplexed `.qza`

### Outputs

- `[run_name]-table.qza`, `[run_name]-repseqs.qza`
- Stats/visualizations under `stats/` including summaries and alpha rarefaction (if enabled)

Proceed to [Taxonomy](taxonomy.md).


