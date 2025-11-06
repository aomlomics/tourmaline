## Configuration

Tourmaline 2 uses three config files, one per step. Example names below reflect defaults; any filename is acceptable.

### 1. QA/QC configuration (config_01_qaqc.yaml)

Key parameters:

```yaml
run_name: your_run
output_dir: /absolute/path/to/results
raw_fastq_path: /abs/path/to/fastqs  # or use sample_manifest_file or trimmed_fastq_path
paired_end: true
to_trim: false

# Trimming (if to_trim)
fwd_primer: ATCG...
rev_primer: ATCG...
discard_untrimmed: false
minimum_length: 100
```

Input options (choose one):

```yaml
raw_fastq_path: /abs/path/to/fastqs
trimmed_fastq_path: /abs/path/to/trimmed
sample_manifest_file: 00-data/manifest_pe.csv  # relative path allowed
```

Manifest formats are documented in [steps/qaqc.md](steps/qaqc.md).

### 2. Repseqs configuration (config_02_repseqs.yaml)

Key parameters:

```yaml
run_name: your_run
output_dir: /absolute/path/to/results
asv_method: dada2pe  # one of: dada2pe, dada2se, deblur

# DADA2 (if dada2*)
dada2_trunc_len_f: 0
dada2pe_trunc_len_r: 0
dada2_trim_left_f: 0
dada2pe_trim_left_r: 0

# Filtering (optional)
to_filter: false
repseq_min_length: 0
repseq_max_length: 100000
repseq_min_abundance: 0.0
repseq_min_prevalence: 0.0
```

Inputs can come from the QA/QC step by matching `run_name`/`output_dir`, or you can provide external `.qza` via:

```yaml
sample_run_name: qa_run_name  # to reuse a different qaqc run name
fastq_qza_file: /abs/path/to/fastq.qza  # external input
```

### 3. Taxonomy configuration (config_03_taxonomy.yaml)

Key parameters:

```yaml
run_name: your_run
output_dir: /absolute/path/to/results
classify_method: naive-bayes  # or consensus-blast, consensus-vsearch, bt2-blca
collapse_taxalevel: 0
classify_threads: 4
```

Repseqs inputs can come from the Repseqs step or be provided externally:

```yaml
repseqs_run_name: repseqs_run  # to reuse a different repseqs run
repseqs_qza_file: /abs/path/to/repseqs.qza
table_qza_file: /abs/path/to/table.qza
```

Reference database parameters:

```yaml
database_name: PR2
refseqs_file: 00-data/refseqs.fna
taxa_file: 00-data/reftax.tsv
pretrained_classifier: /abs/path/to/classifier.qza  # optional for naive-bayes
bowtie_database: /abs/path/to/bt2/index/  # optional for bt2-blca
taxa_ranks: kingdom,phylum,class,order,family,genus,species
```

Method-specific thresholds (examples):

```yaml
# naive-bayes
skl_confidence: 0.7

# consensus methods
perc_identity: 0.8
query_cov: 0.8
min_consensus: 0.51

# bt2-blca
confidence_thres: 0.8
```

See [Running](running.md) for multi-step invocation and [External Data](external_data.md) for conversions.


