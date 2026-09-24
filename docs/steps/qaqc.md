## QA/QC Step

Processes raw FASTQ files (paired or single), provides quality plots, optional primer trimming, optional paired-end merging, and produces QIIME 2 demultiplexed artifacts plus QA summaries.

### Inputs

Choose one in `config_01_qaqc.yaml`:

```yaml
raw_fastq_path: /abs/path/to/raw_fastqs
trimmed_fastq_path: /abs/path/to/trimmed_fastqs
sample_manifest_file: 00-data/manifest_pe.csv
preexisting_fastq_qza: /abs/path/demux.qza
```

Set `paired_end: true|false` accordingly.

### Manifest formats

TSV (current QIIME 2):

Paired-end:

```tsv
sample-id	forward-absolute-filepath	reverse-absolute-filepath
sample1	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz
```

Single-end:

```tsv
sample-id	absolute-filepath
sample1	/path/to/sample1_R1.fastq.gz
```

CSV (legacy):

Paired-end:

```csv
sample-id,absolute-filepath,direction
sample1,/path/to/sample1_R1.fastq.gz,forward
sample1,/path/to/sample1_R2.fastq.gz,reverse
```

Single-end:

```csv
sample-id,absolute-filepath
sample1,/path/to/sample1_R1.fastq.gz
```

### Optional trimming

```yaml
to_trim: true
fwd_primer: ATCG...
rev_primer: ATCG...
discard_untrimmed: false
minimum_length: 100
```

### Optional merging

```yaml
to_merge: true
maxdiffs: 20
merge_stagger: --p-allowmergestagger
```

### Outputs

Written to `[run_name]-qaqc/`:

| File | When | Contents |
|---|---|---|
| `raw_fastq.qza` | `to_trim: true` | Imported reads before trimming |
| `{run_name}_fastq.qza` | always | The artifact downstream steps read |
| `merged_fastq.qza`, `unmerged_fastq.qza` | `to_merge: true` | vsearch merge results |
| `stats/raw_fastq_summary.qzv` | `to_trim: true` | Quality plots before trimming |
| `stats/fastq_summary.qzv` | always | Quality plots of the final reads |
| `stats/merged_fastq_summary.qzv` | `to_merge: true` | Quality plots of merged reads |
| `stats/cutadapt_summary.txt` | `to_trim: true` | Per-sample trimming counts |
| `stats/merge_stats.txt` | `to_merge: true` | Merge success counts |
| `{run_name}-qaqc_config.yaml` | always | Copy of the config that produced this run |

`.qzv` files are QIIME 2 visualizations — open them at [view.qiime2.org](https://view.qiime2.org).

**Read the quality plots before Step 2.** `stats/fastq_summary.qzv` is what you use to choose
DADA2 truncation lengths: find where quality drops off, while leaving enough length for paired
reads to still overlap.

### Optional: per-base quality drop-off report

`check_seq_qual_dropoff` reports, from FastQC/MultiQC output, the position where per-base quality
falls below `seq_quality_cutoff`. It is **not** part of the default `qaqc_all` target — request it
by asking Snakemake for the file:

```bash
snakemake --use-conda -s qaqc_step.Snakefile \
  ../v2-results/my_run-qaqc/stats/my_run-seq_qual_dropoff.txt \
  --configfile config_01_qaqc.yaml --cores 4
```

Continue with [Repseqs](repseqs.md) or run all via [Running](../running.md).


