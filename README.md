<img src="png/tourmaline_banner.png" alt="png/tourmaline_banner" width="100%"/>

[![DOI](https://zenodo.org/badge/125841708.svg)](https://zenodo.org/badge/latestdoi/125841708)

# Tourmaline 2

Tourmaline 2 is an amplicon sequence processing workflow for Illumina sequence data that uses [QIIME 2](https://qiime2.org) and the software packages it wraps. Tourmaline 2 manages commands, inputs, and outputs using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system.

You describe your analysis in YAML configuration files, and Tourmaline runs the right QIIME 2 commands in the right order, keeping a copy of the config next to every set of results so a run can always be traced back to its parameters.

📖 **Full documentation:** [docs/index.md](docs/index.md) (built as an mkdocs site — run `mkdocs serve` locally to browse it)

---

## New to amplicon sequence analysis?

If you are new to metabarcoding or bioinformatics, start here. If you have used QIIME 2 before, skip to [Quick start](#quick-start).

### What this workflow does

You sequenced a marker gene (16S, 18S, COI, 12S MiFish, …) from a set of samples and got back a pile of FASTQ files. Tourmaline turns those into a table of *what organisms were found in which samples*:

```
FASTQ files  →  clean reads  →  unique sequences (ASVs)  →  names for those sequences
   (raw)         Step 1 qaqc      Step 2 repseqs            Step 3 taxonomy
                                       ↓                          ↓
                               feature table + ASV seqs    taxonomy table + barplots
```

### Vocabulary

| Term | Meaning |
|---|---|
| **Amplicon / marker gene** | The short, targeted stretch of DNA you PCR-amplified and sequenced. |
| **Demultiplexed** | One FASTQ file (or pair) per sample, already split out from the pooled sequencing run. Tourmaline expects this. |
| **Primer trimming** | Removing the PCR primer sequences from the start of each read. Done in Step 1 with Cutadapt. |
| **Denoising** | Correcting sequencing errors and collapsing reads into exact unique sequences. Done in Step 2 with DADA2 or Deblur. |
| **ASV** | Amplicon Sequence Variant — one exact unique sequence. The modern replacement for clustered "OTUs". |
| **Representative sequences (repseqs)** | The FASTA-like set of ASV sequences. |
| **Feature table** | A matrix of ASVs (rows) × samples (columns) holding read counts. |
| **Taxonomy assignment** | Comparing each ASV against a reference database to get a name like `Eukaryota;Chordata;…;Lophius americanus`. |
| **Reference database** | Sequences + their known taxonomy (SILVA, PR2, MIDORI, a custom MiFish database, NCBI `nt`, …). |
| **QIIME 2 artifact (`.qza`/`.qzv`)** | A zip file holding data (`.qza`) or a viewable visualization (`.qzv`) plus its provenance. View `.qzv` files at [view.qiime2.org](https://view.qiime2.org). |
| **Conda environment** | An isolated software installation. Tourmaline uses several, by name. |
| **Snakemake rule** | One step of the workflow. Snakemake only re-runs rules whose outputs are missing or out of date. |

### What you need before starting

1. **Demultiplexed FASTQ files** — one file per sample (single-end) or two (paired-end). Tourmaline does not demultiplex.
2. **Your primer sequences** — if you want Tourmaline to trim them.
3. **A reference database** for your marker gene — for taxonomy assignment.
4. *(Optional)* **A sample metadata TSV** — used for barplots and diversity plots.

### Practical advice for a first run

- Start with the small example data in [`00-data/`](00-data/) to confirm your install works before pointing at real data.
- Run one step at a time and look at the `.qzv` outputs before moving on. In particular, look at the quality plots from Step 1 *before* choosing DADA2 truncation lengths in Step 2.
- Give each attempt a distinct `run_name`. Tourmaline keeps runs side by side, so comparing parameter sets is cheap.
- Use `--dryrun` (see [Running](#running-the-workflow)) to see what *would* run without running it.

---

## Quick start

```bash
# 1. Get Tourmaline (the default branch, V2)
git clone https://github.com/aomlomics/tourmaline.git
cd tourmaline

# 2. Create the environments (one time) — see Setup below for the full list
conda create -c conda-forge -c bioconda -n snakemake-tour2 snakemake biopython yq parallel

# 3. Activate the Snakemake environment
conda activate snakemake-tour2

# 4. Edit the example configs, then run all three steps
./tourmaline.sh \
  --step qaqc,repseqs,taxonomy \
  --configfile config_01_qaqc.yaml,config_02_repseqs.yaml,config_03_taxonomy.yaml \
  --cores 6
```

See [docs/quick_start.md](docs/quick_start.md) and [docs/install.md](docs/install.md) for more.

---

## Overview of the steps

Tourmaline 2 is modular. Each step has its own Snakefile and its own config file, and steps chain together through files on disk — so **any step can be the starting point** if you supply correctly formatted input.

### Step 1 — Sequence QA/QC (`qaqc`)

* Processes demultiplexed FASTQ files (paired-end or single-end).
* Optionally trims primer sequences from raw reads (Cutadapt).
* Optionally merges paired-end reads with vsearch (`to_merge`), e.g. when you plan to use Deblur.
* Produces sequence quality plots for raw and/or trimmed reads.
* Creates a QIIME 2 demultiplexed sequence artifact.

📖 [docs/steps/qaqc.md](docs/steps/qaqc.md)

### Step 2 — Representative sequences (`repseqs`)

* Generates ASVs with DADA2 (paired- or single-end) or Deblur (single-end).
* Optional filtering by length, abundance, prevalence, frequency, and sample count.
* Produces the feature table and representative sequences, plus denoising stats.
* Optional alpha rarefaction and core diversity metrics (`plot_diversity`).

📖 [docs/steps/repseqs.md](docs/steps/repseqs.md)

### Step 3 — Taxonomy assignment (`taxonomy`)

* Assigns taxonomy using one of five methods:
  * [Naive Bayes classifier as implemented in QIIME 2](https://docs.qiime2.org/2024.10/plugins/available/feature-classifier/classify-sklearn/)
  * [Consensus BLAST as implemented in QIIME 2](https://docs.qiime2.org/2024.10/plugins/available/feature-classifier/classify-consensus-blast/)
  * [Consensus VSEARCH as implemented in QIIME 2](https://docs.qiime2.org/2024.10/plugins/available/feature-classifier/classify-consensus-vsearch/)
  * [Anacapa's Bowtie 2 and BLCA method](https://github.com/limey-bean/Anacapa?tab=readme-ov-file#step-3-taxonomic-assignment-using-bowtie-2-and-blca)
  * [REVAMP's BLASTn against NCBI nt with lowest common ancestor](https://github.com/McAllister-NOAA/REVAMP)
* Produces a taxonomy table, an interactive taxa barplot, a table collapsed to a chosen rank, and a combined ASV/taxonomy/sequence TSV.
* Optional interactive [Krona](https://github.com/marbl/Krona/wiki) plot (`make_krona`), for any classify method.

📖 [docs/steps/taxonomy.md](docs/steps/taxonomy.md)

### Analysis metadata (utility script)

* Creates a file with metadata about the analysis using FAIR eDNA terms.
* File can be read into the [NOAA Ocean DNA Explorer](https://www.ngi.msstate.edu/node).
* This is a standalone script, not a `tourmaline.sh` step — see [Generating analysis metadata](#generating-analysis-metadata).

📖 [docs/metadata.md](docs/metadata.md)

---

## Setup

### Required

* [Conda (Miniconda works well)](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
* [QIIME 2 (2024.10) amplicon distribution](https://docs.qiime2.org/2024.10/install/) — **the environment must be named exactly `qiime2-amplicon-2024.10`**, because the Snakemake rules request it by name.
* A Snakemake environment with a few extra packages:

   ```bash
   conda create -c conda-forge -c bioconda -n snakemake-tour2 snakemake biopython yq parallel
   ```

* Tourmaline itself (the default branch is `V2`):

   ```bash
   git clone https://github.com/aomlomics/tourmaline.git
   ```

### Optional, per feature

| Environment | Needed for | Create with |
|---|---|---|
| `bt2-blca` | `classify_method: bt2-blca` | `conda create -c conda-forge -c bioconda -n bt2-blca biopython muscle=3.8 bowtie2` |
| `revamp` | `classify_method: revamp` | `conda create -c conda-forge -c bioconda -n revamp "blast>=2.13" "taxonkit>=0.20" r-base r-dplyr bioconductor-biostrings perl perl-list-moreutils krona` |
| `krona` | `make_krona: True` | `conda create -c conda-forge -c bioconda -n krona krona` |

REVAMP additionally needs a clone of [REVAMP](https://github.com/McAllister-NOAA/REVAMP) and a local NCBI `nt` BLAST database with prepared taxonomy files — see [docs/steps/taxonomy.md#revamp](docs/steps/taxonomy.md#revamp).

📖 [docs/install.md](docs/install.md)

### Running requirements

* The `snakemake-tour2` environment must be **activated**.
* You must run from the Tourmaline directory (it contains `tourmaline.sh`, the Snakefiles, and `scripts/`, which rules call by relative path). Output can go anywhere via `output_dir`.
* Config files for each step you intend to run.

---

## Configuration files

Each step reads one YAML config file. The example configs in the repository are extensively commented and are the best starting point — copy one and edit it:

| Step | Example config | Snakefile |
|---|---|---|
| qaqc | [`config_01_qaqc.yaml`](config_01_qaqc.yaml) | `qaqc_step.Snakefile` |
| repseqs | [`config_02_repseqs.yaml`](config_02_repseqs.yaml) | `repseqs_step.Snakefile` |
| taxonomy | [`config_03_taxonomy.yaml`](config_03_taxonomy.yaml) | `taxonomy_step.Snakefile` |

Config files can have any name; pass whichever you want with `--configfile`.

> **Keep optional keys present but empty.** The Snakefiles read most config values directly, so a *missing* key raises a `KeyError` before the workflow starts. An empty value is fine and means "not set" — deleting the line is not.

📖 **[docs/configuration.md](docs/configuration.md) is the complete parameter reference.** The summaries below cover the keys most runs need.

### 1. QA/QC configuration

```yaml
run_name: my_run                # prefix for this run's outputs
output_dir: "../v2-results"     # where outputs are written
paired_end: True                # True for paired-end, False for single-end
to_trim: True                   # trim primers with Cutadapt
to_merge: False                 # merge paired reads with vsearch
assay_name: Bacteria-16S-V4V5-Parada   # for metadata reporting

# Primer trimming (used when to_trim: True)
fwd_primer: GTGYCAGCMGCCGCGGTAA # IUPAC ambiguity codes supported
rev_primer: GGACTACNVGGGTWTCTAAT
discard_untrimmed: False        # drop reads with no primer match
minimum_length: 50              # minimum read length kept after trimming
trimming_threads: 5

# Merging (used when to_merge: True)
maxdiffs: 20
merge_stagger: --p-allowmergestagger
```

`assay_name` should come from the [NOAA Omics metabarcoding assays](https://github.com/NOAA-Omics/noaa-omics-metabarcoding-assays/blob/main/assays.tsv) controlled vocabulary.

#### QA/QC input files

Choose **one** of these four and leave the others blank:

```yaml
# Directory of raw demultiplexed fastq files. Sample names are the file name prefixes.
raw_fastq_path: [path]
# Directory of already-trimmed fastq files (suffix _R[1,2].fastq.gz).
trimmed_fastq_path: [path]
# A QIIME 2 manifest file. Can point to trimmed or untrimmed reads.
sample_manifest_file: [path/filename]
# An already-imported QIIME 2 demultiplexed sequence artifact.
preexisting_fastq_qza: [path]
```

#### FASTQ file naming (when not using a manifest)

* Paired-end: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
* Alternative: `{sample}_R1_001.fastq.gz` and `{sample}_R2_001.fastq.gz`
* Single-end: `{sample}_R1.fastq.gz` or `{sample}_R1_001.fastq.gz`

#### Sample manifest format

Provide either the current QIIME 2 tab-separated format or the legacy comma-separated format. Headers must match exactly.

**Tab-separated**, paired-end:

```tsv
sample-id	forward-absolute-filepath	reverse-absolute-filepath
sample1	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz
```

**Tab-separated**, single-end:

```tsv
sample-id	absolute-filepath
sample1	/path/to/sample1_R1.fastq.gz
```

**CSV (legacy)**, paired-end:

```csv
sample-id,absolute-filepath,direction
sample1,/path/to/sample1_R1.fastq.gz,forward
sample1,/path/to/sample1_R2.fastq.gz,reverse
```

**CSV (legacy)**, single-end:

```csv
sample-id,absolute-filepath
sample1,/path/to/sample1_R1.fastq.gz
```

### 2. Representative sequences configuration

```yaml
run_name: my_run                # can match or differ from the qaqc run
output_dir: "../v2-results"
asv_method: dada2pe             # dada2pe | dada2se | deblur
asv_threads: 5

# DADA2 parameters (dada2pe / dada2se)
dada2_trunc_len_f: 245          # truncate forward reads here (0 = full length)
dada2pe_trunc_len_r: 190        # truncate reverse reads (paired-end only)
dada2_trim_left_f: 0            # bases trimmed from the start of forward reads
dada2pe_trim_left_r: 0          # bases trimmed from the start of reverse reads
dada2_max_ee_f: 2               # max expected errors, forward
dada2pe_max_ee_r: 2             # max expected errors, reverse
dada2_trunc_q: 2                # truncate at the first base with this quality
dada2_pooling_method: pseudo    # independent | pseudo | pooled
dada2_chimera_method: consensus # none | consensus | pooled

# Deblur parameters (deblur)
deblur_trim_length: 150         # final sequence length
reference_seqs:                 # reference artifact for positive filtering

# Optional diversity outputs
plot_diversity: True
alpha_max_depth: 500            # required when plot_diversity: True
core_sampling_depth: 500        # required when plot_diversity: True

# Optional filtering
to_filter: False
repseq_min_length: 0
repseq_max_length: 0
repseq_min_abundance: 0         # relative abundance, 0–1
repseq_min_prevalence: 0        # fraction of samples, 0–1
repseq_min_frequency: 0         # absolute total count
repseq_min_samples: 0           # minimum number of samples
```

> ⚠️ **If you set `to_filter: True`, also set `repseq_max_length` to a real upper bound** (e.g. `100000`). The length filter keeps sequences where `length <= repseq_max_length`, so leaving the default `0` removes every sequence.

Choosing truncation lengths is the main judgement call in this step: look at the quality plots produced by Step 1, and truncate where quality drops off while still leaving enough overlap for paired reads to merge.

#### Repseqs input files

Two options:

**1) Use an existing Tourmaline QA/QC run**

* Use the same `run_name` and `output_dir` for both steps, or
* Use a different `run_name` for repseqs and set `qaqc_run_name` to the QA/QC run you want. Useful when testing different trimming parameters against one denoising setup.

**2) Provide an externally generated QIIME 2 sequence artifact**

Set `fastq_qza_file`. See [Starting from external data](#starting-from-external-data).

### 3. Taxonomy configuration

```yaml
run_name: my_run
output_dir: "../v2-results"
classify_method: naive-bayes    # naive-bayes | consensus-blast | consensus-vsearch | bt2-blca | revamp
collapse_taxalevel: 7           # rank (1–7) for the additional collapsed count table
classify_threads: 10
sample_metadata_file: 00-data/metadata.tsv   # optional, for the barplot
```

Reference database keys:

```yaml
database_name: "silva-138_1-99-515f_926r-uniq"  # descriptive only, recorded in metadata
refseqs_file: [path]            # reference sequences (.qza or FASTA)
taxa_file: [path]               # reference taxonomy (.qza or TSV)
taxa_ranks: kingdom,phylum,class,order,family,genus,species  # must match the database
pretrained_classifier: [path]   # naive-bayes only; overrides refseqs_file/taxa_file
bowtie_database: [path]         # bt2-blca only; prebuilt index, else built from refseqs
```

Method-specific parameters:

```yaml
# naive-bayes
skl_confidence: 0.7      # confidence threshold limiting assignment depth

# consensus-blast / consensus-vsearch
perc_identity: 0.8       # minimum percent identity for a hit (0–1)
query_cov: 0.8           # minimum query coverage for a hit (0–1)
min_consensus: 0.51      # fraction of hits that must agree

# bt2-blca
confidence_thres: 0.8    # bootstrap confidence threshold limiting assignment depth

# revamp (see docs/steps/taxonomy.md#revamp for database setup)
revamp_dir: [path]              # clone of https://github.com/McAllister-NOAA/REVAMP
revamp_blastdb: [path]          # NCBI nt directory with a taxdump/ from ncbi_db_cleanup.sh
revamp_blast_results: [path]    # optional: a BLASTn btab produced elsewhere
revamp_blast_mode: mostEnvOUT   # allIN | allEnvOUT | mostEnvOUT
revamp_query_cov: 90            # percent of ASV length a hit must cover (0–100)
revamp_taxonomy_cutoffs: "97,95,90,80,70,60"  # percent ID cutoffs, ordered S,G,F,O,C,P

# Krona plot (any classify method)
make_krona: False        # writes figures/{run_name}-krona.html
krona_per_sample: True   # add one Krona dataset per sample

# Extra flags appended to the classifier command
classify_params: --verbose
```

#### Taxonomy input files

Two options:

**1) Use an existing Tourmaline repseqs run**

* Use the same `run_name` and `output_dir` for both steps, or
* Use a different `run_name` for taxonomy and set `repseqs_run_name`. Useful when comparing classifiers against one set of ASVs.

**2) Provide externally generated QIIME 2 artifacts**

Set both `repseqs_qza_file` and `table_qza_file`. See [Starting from external data](#starting-from-external-data).

---

## Starting from external data

Unlike Tourmaline 1, you can begin at any step using data from another program, as long as it is formatted as the QIIME 2 artifact that step expects. For example, if you already have ASV sequences and only want taxonomy, import them and set `repseqs_qza_file`.

Activate the QIIME 2 environment first:

```bash
conda activate qiime2-amplicon-2024.10
```

**Demultiplexed reads → `.qza`** (for `preexisting_fastq_qza` or `fastq_qza_file`). Needs a [manifest file](#sample-manifest-format):

```bash
# Paired-end
qiime tools import \
   --type 'SampleData[PairedEndSequencesWithQuality]' \
   --input-path my_pe.manifest \
   --output-path output-file_pe_fastq.qza \
   --input-format PairedEndFastqManifestPhred33V2

# Single-end
qiime tools import \
   --type 'SampleData[SequencesWithQuality]' \
   --input-path my_se.manifest \
   --output-path output-file_se_fastq.qza \
   --input-format SingleEndFastqManifestPhred33V2
```

**ASV sequences (FASTA) → `.qza`** (for `repseqs_qza_file`):

```bash
qiime tools import \
   --type 'FeatureData[Sequence]' \
   --input-path my-asvs.fasta \
   --output-path output-asvs.qza
```

**Read count table → `.qza`** (for `table_qza_file`). If you have a BIOM file, [check its format first](https://docs.qiime2.org/2024.10/tutorials/importing/#feature-table-data):

```bash
# BIOM v1.0.0
qiime tools import \
  --input-path feature-table-v100.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path feature-table.qza
```

If you have a TSV with unique sequences as rows and samples as columns, [convert to BIOM](https://biom-format.org/documentation/biom_conversion.html) first:

```bash
biom convert -i otu_table.txt -o new_otu_table.biom --to-hdf5 --table-type="OTU table"

qiime tools import \
  --input-path new_otu_table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path feature-table.qza
```

> QIIME 2 artifacts written by one QIIME 2 version are often **not** readable by earlier versions. Import with 2024.10 if you plan to run Tourmaline 2 on the result.

📖 [docs/external_data.md](docs/external_data.md)

---

## Running the workflow

### Basic usage

Activate the environment and run from the Tourmaline directory:

```bash
conda activate snakemake-tour2
./tourmaline.sh --step/-s [steps] --configfile/-c [config_files] --cores/-n [num_cores]
```

Run a single step:

```bash
./tourmaline.sh -s taxonomy -c config_03_taxonomy.yaml -n 6
```

Run all steps with one command:

```bash
./tourmaline.sh -s qaqc,repseqs,taxonomy \
  -c config_01_qaqc.yaml,config_02_repseqs.yaml,config_03_taxonomy.yaml -n 6
```

**Important:**

* The number of steps must match the number of config files.
* Config files must be given in the same order as the steps.
* Valid steps are `qaqc`, `repseqs`, and `taxonomy`.

### Running Snakemake directly

You can still call Snakemake yourself — necessary for dry runs, single rules, or `--printshellcmds`. Each step has its own Snakefile and target rule:

```bash
snakemake --use-conda -s qaqc_step.Snakefile     qaqc_all     --configfile config_01_qaqc.yaml --cores 6 --dryrun
snakemake --use-conda -s repseqs_step.Snakefile  run_denoise  --configfile config_02_repseqs.yaml --cores 6
snakemake --use-conda -s taxonomy_step.Snakefile run_taxonomy --configfile config_03_taxonomy.yaml --cores 6
```

A dry run (`--dryrun`) is the fastest way to check a config before committing compute to it.

### Parameter sweeps

To compare many parameter sets, expand a base config over a parameter space and run the results in parallel:

```bash
python scripts/generate_configs.py <base_config> <parameter_space_config>

scripts/run_parallel_tourmaline.sh \
  --config-dir parameter_sweep_configs --config-prefix config-01-qaqc \
  --step qaqc --parallel-jobs 4 --cores-per-job 6
```

Example parameter space files are in [`00-data/`](00-data/) (`parameter_space*.yaml`). Parallel runs require GNU `parallel`.

### HPC / SLURM

`scripts/sbatch_tourmaline2_step*.sh` are SLURM array wrappers for the same steps. Adapt the account, partition, and resource lines to your cluster.

📖 [docs/running.md](docs/running.md)

---

## Outputs

Everything lands under `output_dir`, one directory per run and step:

```
output_dir/
├── [run_name]-qaqc/       # QA/QC outputs
├── [run_name]-repseqs/    # Representative sequences outputs
└── [run_name]-taxonomy/   # Taxonomy assignment outputs
```

Each step also **copies its config file into its own output directory** as `{run_name}-{step}_config.yaml`, so a run's provenance sits next to its results.

Key files:

| File | Step | Contents |
|---|---|---|
| `raw_fastq.qza`, `{run_name}_fastq.qza` | qaqc | Imported reads, before and after trimming |
| `stats/*.qzv` | qaqc | Quality summaries — view at [view.qiime2.org](https://view.qiime2.org) |
| `{run_name}-table.qza` / `-table.tsv` | repseqs | Feature table (ASV × sample counts) |
| `{run_name}-repseqs.qza` | repseqs | Representative sequences |
| `stats/dada2_stats.tsv` / `deblur_stats.tsv` | repseqs | Reads retained at each denoising stage |
| `{run_name}-taxonomy.tsv` / `.qza` | taxonomy | Per-ASV taxonomy assignments |
| `figures/{run_name}-taxa_barplot.qzv` | taxonomy | Interactive taxa barplot |
| `{run_name}-taxa_sample_table_l{N}.tsv` | taxonomy | Counts collapsed to rank N |
| `{run_name}-asv_taxa_features.tsv` | taxonomy | Combined ASV + taxonomy + sequence table |
| `figures/{run_name}-krona.html` | taxonomy | Krona plot, when `make_krona: True` |

---

## Generating analysis metadata

To produce a FAIR eDNA analysis metadata TSV (readable by the [NOAA Ocean DNA Explorer](https://www.ngi.msstate.edu/node)), run `scripts/format_analysisMetadata.py` **after** your runs finish. It takes the results directory and the **run names** of each step — it finds each run's copied config file itself — plus a `project_id` and an output folder.

```bash
python scripts/format_analysisMetadata.py \
  -w ../v2-results \
  -q my_qaqc_run -r my_repseqs_run -t my_taxonomy_run \
  -p my_project \
  -O output_folder/
```

```
usage: format_analysisMetadata.py [-h] -w WORKING_DIR -q QAQC_RUN_NAME -r REPSEQS_RUN_NAME
                                  -t TAXONOMY_RUN_NAME -p PROJECT_ID [-a ASSAY_NAME]
                                  [-A ANALYSIS_RUN_NAME] [-T TOURMALINE_METADATA] -O OUTPUT_FOLDER

options:
  -h, --help            show this help message and exit
  -w, --working_dir     Working directory containing the step output folders
  -q, --qaqc_run_name   Run name for the qaqc step
  -r, --repseqs_run_name
                        Run name for the repseqs step
  -t, --taxonomy_run_name
                        Run name for the taxonomy step
  -p, --project_id      Value for project_id
  -a, --assay_name      Value for assay_name, otherwise uses value in qaqc config
  -A, --analysis_run_name
                        Value for analysis_run_name, otherwise uses taxonomy run name
  -T, --tourmaline_metadata
                        Path to tourmaline metadata (default ./00-data/tourmaline_metadata.yaml)
  -O, --output_folder   Output folder where files will be saved
```

Alongside the metadata TSV it copies the taxonomy and table outputs into the output folder with the analysis run name as a prefix.

📖 [docs/metadata.md](docs/metadata.md)

---

## Documentation map

| Page | Covers |
|---|---|
| [docs/index.md](docs/index.md) | Overview and feature summary |
| [docs/quick_start.md](docs/quick_start.md) | Shortest path to a first run |
| [docs/install.md](docs/install.md) | Requirements and conda environments |
| [docs/configuration.md](docs/configuration.md) | **Complete** parameter reference for every config file |
| [docs/running.md](docs/running.md) | `tourmaline.sh`, direct Snakemake, sweeps, HPC |
| [docs/steps/qaqc.md](docs/steps/qaqc.md) | Step 1 details |
| [docs/steps/repseqs.md](docs/steps/repseqs.md) | Step 2 details |
| [docs/steps/taxonomy.md](docs/steps/taxonomy.md) | Step 3 details, including REVAMP and Krona |
| [docs/external_data.md](docs/external_data.md) | Starting from externally generated inputs |
| [docs/metadata.md](docs/metadata.md) | FAIR eDNA analysis metadata |
| [docs/troubleshooting.md](docs/troubleshooting.md) | Common errors and fixes |
| [docs/citation_legacy.md](docs/citation_legacy.md) | Citation and v1 resources |

---

## Major changes in v2 vs. v1

**To use the legacy v1 version of Tourmaline**, check out the [V1 branch](https://github.com/aomlomics/tourmaline/tree/V1) of this repository. The v1 README is kept here as [`README_v1.md`](README_v1.md).

* **Run via the `tourmaline.sh` script.** Instead of interacting with Snakemake rules directly, the main entry point is `tourmaline.sh`, which runs one or more steps, each with its own config file and a shared core count. You can still call individual Snakemake rules — each step has its own Snakefile, so specify the right one.
* **Modular steps.** qaqc, repseqs, and taxonomy are separate Snakefiles with separate configs and separate output directories, so you can re-run one step over several parameter sets while reusing upstream output.
* **Externally-generated data.** Any step can start from data produced elsewhere, as long as it is formatted as the expected QIIME 2 artifact.
* **More taxonomy methods.** bt2-blca and REVAMP join the three QIIME 2 classifiers.
* **FAIR eDNA metadata output** for submission to the NOAA Ocean DNA Explorer.

📖 [docs/citation_legacy.md](docs/citation_legacy.md)

## Citation

Thompson, L. R., Anderson, S. R., Den Uyl, P. A., Patin, N. V., Sanderson, G. & Goodwin, K. D. Tourmaline: A containerized workflow for rapid and iterable amplicon sequence analysis using QIIME 2 and Snakemake. *GigaScience*, Volume 11, 2022, giac066. doi: [10.1093/gigascience/giac066](https://doi.org/10.1093/gigascience/giac066)

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
