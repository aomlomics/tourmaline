## Configuration

Tourmaline 2 uses three config files, one per step. Example names below reflect defaults; any filename is acceptable. Templates for each file are in the Tourmaline folder.

### 1. QA/QC configuration (Template: config_01_qaqc.yaml)

The qaqc config defines how Tourmaline imports, trims, and summarizes raw reads.

**Core run settings (required)**

```yaml
run_name: my_run                 # unique identifier; used in output folders
output_dir: /path/to/results     # absolute or relative directory for outputs
paired_end: true                 # true for paired-end, false for single-end data
to_trim: true                    # enable primer trimming with Cutadapt
to_merge: false                  # enable vsearch merge (paired-end only)
to_filter: false                # enable feature filtering steps
assay_name: Bacteria-16S-V4V5-Parada  # for metadata reporting
```
Select assay name based on the [NOAA Omics metabarcoding assays](https://github.com/NOAA-Omics/noaa-omics-metabarcoding-assays/blob/main/assays.tsv) controlled vocabulary. If your assay is not available, please create an [issue](https://github.com/NOAA-Omics/noaa-omics-metabarcoding-assays/issues).

**Input data sources (choose one)**

```yaml
raw_fastq_path: /abs/path/to/raw_fastqs       # directory of raw FASTQ(.gz)
trimmed_fastq_path: /abs/path/to/trimmed      # directory of already trimmed FASTQs
sample_manifest_file: 00-data/manifest.tsv    # QIIME 2 manifest (TSV or CSV)
preexisting_fastq_qza: 00-data/demux.qza      # existing demuxed artifact
```

Provide only the fields relevant to your data source; leave others blank. Manifest formats are detailed in [steps/qaqc.md](steps/qaqc.md).

**Primer trimming parameters (required when `to_trim: true`)**

```yaml
fwd_primer: GTGYCAGCMGCCGCGGTAA      # IUPAC supported
rev_primer: GGACTACNVGGGTWTCTAAT
discard_untrimmed: false             # discard reads without primer match
minimum_length: 50                   # post-trimming minimum length (bp)
trimming_threads: 5                  # threads for Cutadapt
```

**Merging + compute options (required when `to_merge: true`)**

```yaml
maxdiffs: 20                         # vsearch merge mismatches
merge_stagger: --p-allowmergestagger # optional vsearch flag
```

### 2. Repseqs configuration (Template: config_02_repseqs.yaml)

Controls ASV generation, filtering, and optional diversity plots.

**Core run settings (required)**

```yaml
run_name: my_run
output_dir: /path/to/results
asv_method: dada2pe          # one of: dada2pe | dada2se | deblur
asv_threads: 5               # threads passed to denoisers
```

**Input data sources (choose one, otherwise will default to** \[output_dir\]/\[my_run-qaqc\])

```yaml
qaqc_run_name: my_qaqc_run         # reuse QA/QC outputs from another run
fastq_qza_file: /abs/path/demux.qza  # external demultiplexed sequences
```

If neither is supplied, the workflow expects demultiplexed reads from the QA/QC step with the same `run_name`.

**Metadata + diversity options**

```yaml
sample_metadata_file: 00-data/metadata.tsv   # optional metadata for summaries/diversity
plot_diversity: true                         # produce alpha/core metrics outputs
alpha_max_depth: 500                         # required when plot_diversity is true
core_sampling_depth: 500                     # required when plot_diversity is true
```

**DADA2 parameters (required when `asv_method` starts with dada2)**

```yaml
dada2_trunc_len_f: 245        # forward truncation length
dada2pe_trunc_len_r: 190      # reverse truncation (paired-end only)
dada2_trim_left_f: 0          # forward trim from left
dada2pe_trim_left_r: 0        # reverse trim from left
dada2_max_ee_f: 2             # forward max expected errors
dada2pe_max_ee_r: 2           # reverse max expected errors
dada2_trunc_q: 2              # truncate at quality score
dada2_pooling_method: pseudo  # independent | pseudo | pooled
dada2_chimera_method: consensus
dada2_min_fold_parent_over_abundance: 1
dada2_n_reads_learn: 1000000
dada2_hashed_feature_ids: --p-hashed-feature-ids  # optional
```

**Deblur parameters (required when `asv_method: deblur`)**

```yaml
deblur_trim_length: 150       # final sequence length (bp)
deblur_trim_left: 0
deblur_mean_error: 0.005
deblur_min_reads: 2
deblur_min_size: 2
deblur_indel_max: 3
reference_seqs: 00-data/ref.qza  # required reference set
```

**Post-denoising filtering (required if `to_filter` is `True`)**

```yaml
repseq_min_length: 0
repseq_max_length: 0
repseq_min_abundance: 0
repseq_min_prevalence: 0
repseq_min_frequency: 0
repseq_min_samples: 0
```

### 3. Taxonomy configuration (Template: config_03_taxonomy.yaml)

Defines how representative sequences are assigned taxonomy and summarized.

**Core run settings (required)**

```yaml
run_name: my_run
output_dir: /path/to/results
classify_method: naive-bayes      # options: naive-bayes | consensus-blast | consensus-vsearch | bt2-blca
taxa_ranks: kingdom,phylum,class,order,family,genus,species
collapse_taxalevel: 7             # taxonomy level for collapsed table
classify_threads: 10
```

**Input data sources (choose one)**

```yaml
repseqs_run_name: my_repseqs_run          # reuse outputs from another run
repseqs_qza_file: /abs/path/repseqs.qza   # external representative sequences
table_qza_file: /abs/path/table.qza       # external feature table
```

If no external inputs are supplied, the workflow uses artifacts produced by the Repseqs step with matching `run_name`.

**Reference database parameters**

```yaml
database_name: silva-138_1
refseqs_file: 00-data/silva-seqs.qza    # required unless using pretrained classifier
taxa_file: 00-data/silva-tax.qza        # required unless using pretrained classifier
sample_metadata_file: 00-data/metadata.tsv  # optional for barplots
```

**Naive Bayes options**

```yaml
pretrained_classifier: /abs/path/classifier.qza  # optional, overrides refseqs/taxa files
skl_confidence: 0.7                              # confidence threshold
```

**Consensus BLAST/VSEARCH options**

```yaml
perc_identity: 0.8
query_cov: 0.8
min_consensus: 0.51
```

**BT2-BLCA options**

```yaml
bowtie_database: /abs/path/bowtie2_index/   # optional; auto-built if omitted
confidence_thres: 0.8
```

**Additional classifier flags**

```yaml
classify_params: --verbose   # appended to the chosen classifier command
```

See [Running](running.md) for multi-step invocation and [External Data](external_data.md) for conversions and artifact preparation tips.

### 4. Tax-credit configuration (Template: config_04_tax_credit.yaml)

The tax-credit step benchmarks reference databases using [tax-credit](../tax-credit/) simulations and Tourmaline taxonomy assignment rules. It does not require outputs from the repseqs or taxonomy steps.

**Core run settings**

```yaml
run_name: mifish_tax_credit
output_dir: ../v2-results
tax_credit_package_dir: ../tax-credit   # pip install -e this path in qiime2 env
```

**Reference databases (one or more)**

```yaml
reference_databases:
  - id: my_database
    refseqs_file: /path/to/sequences.fasta   # or .qza
    taxa_file: /path/to/taxonomy.txt         # or .qza
    fwd_primer: GCCGGTAAAACTCGTGCCAGC
    rev_primer: CATAGTGGGGTATCTAATCCCAGTTTG
    fwd_primer_id: MiFishF
    rev_primer_id: MiFishR
    read_length: 250            # required when truncate is true
    min_read_length: 140
    trim_primers: true          # false for pre-trimmed amplicon references (e.g. rCRUX)
    truncate: true
```

**Evaluation methods (one or more)**

```yaml
evaluation_methods:
  - cross-validated          # taxonomy-aware CV folds
  - cross-validated-trad     # traditional KFold CV
  - novel-taxa               # novel-taxa simulation
  - self-validated           # full database classified against itself
  # - mock-community         # see mock_communities below
```

**Simulation parameters**

```yaml
iterations: 10
novel_taxa_levels: [6, 5, 4, 3]
cv_recall_max_level: 6
novel_recall_min_level: 3
force_regenerate: false
```

`cv_recall_max_level` controls cross-validated assignment manifest generation (default `6`; `min_level` is always `max_level - 1` so each CV fold is listed once). Per-database simulation settings (`read_length`, `min_read_length`, `trim_primers`, `truncate`) are configured on each `reference_databases` entry (see above).

**Taxonomic assignment** — uses the same keys as the taxonomy step (`classify_method`, `skl_confidence`, `classify_params`, etc.). Assignment runs via Snakemake rules shared with `taxonomy_step.Snakefile`, not tax-credit shell templates.

```yaml
classify_method: naive-bayes
classify_threads: 5
nb_confidence_values: [0.7]
blca_confidence_values: [0.8]
skl_confidence: 0.7
confidence_thres: 0.8
fit_params: "--p-feat-ext--ngram-range '[7,7]' --p-classify--alpha 0.001"
generate_plots: true
```

Assignment jobs write method-relevant parameters to `assignment_manifest.tsv`; unused fields are left blank. List-valued `perc_identity`, `query_cov`, and `min_consensus` expand into parameter sweeps for consensus methods (and `perc_identity` / `query_cov` for bt2-blca).

**Mock community** (when `mock-community` is listed in `evaluation_methods`)

```yaml
mock_communities:
  - id: fuhrman_18Sv4
    feature_table_biom: /path/to/feature_table.biom
    rep_seqs_fasta: /path/to/rep_seqs.fna
    references:
      - id: pr2-ssu
        expected_dir: /path/to/expected
```

Install tax-credit in the QIIME 2 amplicon environment before running:

```bash
pip install -e ../tax-credit
```


