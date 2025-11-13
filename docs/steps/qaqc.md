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

- Demultiplexed sequences (`raw_fastq.qza` or trimmed) and summary visualizations
- Stats under `[run_name]-qaqc/stats/` (raw summaries, trimmed summaries, merge stats when enabled, sequence-quality drop-off report)

Continue with [Repseqs](repseqs.md) or run all via [Running](../running.md).


