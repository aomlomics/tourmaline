## Tourmaline 2 Documentation

Tourmaline 2 is an amplicon sequence processing workflow for Illumina sequence data that uses [QIIME 2](https://qiime2.org) and the software packages it wraps. Tourmaline 2 manages commands, inputs, and outputs using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system.

**Tourmaline 2 uses [QIIME 2 2024.10](https://docs.qiime2.org/2024.10/install/) amplicon workflow.** To use the Legacy V1 version of Tourmaline, check out the [V1 branch](https://github.com/aomlomics/tourmaline/tree/V1) of this repository.

## Why should I use Tourmaline?

Tourmaline has several features that enhance usability and interoperability:

* **Portability.** Native support for Linux and macOS in addition to Docker containers.
* **QIIME 2.** The core commands of Tourmaline, including the [DADA2](https://benjjneb.github.io/dada2/index.html) and [Deblur](https://github.com/biocore/deblur) packages, are all commands of QIIME 2, one of the most popular amplicon sequence analysis software tools available. You can print all of the QIIME 2 and other shell commands of your workflow before or while running the workflow.
* **Snakemake.** Managing the workflow with Snakemake provides several benefits: 
  - **Configuration file** contains all parameters in one file, so you can see what your workflow is doing and make changes for a subsequent run.
  - **Directory structure** is the same for every Tourmaline run, so you always know where your outputs are.
  - **On-demand commands** mean that only the commands required for output files not yet generated are run, saving time and computation when re-running part of a workflow.
* **Parameter optimization.** The configuration file and standard directory structure make it simple to test and compare different parameter sets to optimize your workflow. Included code helps choose read truncation parameters and identify outliers in representative sequences (ASVs).
* **Visualizations and reports.** Every Tourmaline run produces visualizations and summaries with links to web-viewable QIIME 2 visualization files.
* **Downstream analysis.** Analyze the output of single or multiple Tourmaline runs programmatically, with qiime2R in R or the QIIME 2 Artifact API in Python, using the provided R and Python notebooks or your own code.

## What QIIME 2 options does Tourmaline support?

If you have used QIIME 2 before, you might be wondering which QIIME 2 commands Tourmaline uses and supports. All commands are specified as rules in the Snakefiles, and typical workflows are shown as directed acyclic graphs. The main analysis features and options supported by Tourmaline are as follows:

* FASTQ sequence import using a manifest file, a folder of fastq.gz files, or use your pre-imported FASTQ .qza file
* Denoising with [DADA2](https://doi.org/10.1038/nmeth.3869) (paired-end and single-end) and [Deblur](https://doi.org/10.1128/msystems.00191-16) (single-end)
* Feature classification (taxonomic assignment) with options of [naive Bayes](https://doi.org/10.1186/s40168-018-0470-z), consensus [BLAST](https://doi.org/10.1186/1471-2105-10-421), consensus [VSEARCH](https://doi.org/10.7717/peerj.2584), and [BT2-BLCA](https://github.com/limey-bean/Anacapa?tab=readme-ov-file#step-3-taxonomic-assignment-using-bowtie-2-and-blca)
* Feature filtering by taxonomy, sequence length, feature ID, and abundance/prevalence
* Interactive taxonomy barplots and visualizations
* Alpha diversity metrics, rarefaction analyses, and ordination plots

## Major changes in v2 vs. v1

### Run via tourmaline.sh script

Instead of interacting with Snakemake rules directly, the main way to run Tourmaline 2 is through the `tourmaline.sh` script. This script allows you to run one or more of the workflow steps at a time, specify specific config files, and set the maximum number of cores. You must be located in the tourmaline directory when running it, however you can set the output file destinations to anywhere.

Usage:

```bash
conda activate snakemake-tour2
./tourmaline.sh --step [qaqc,repseqs,taxonomy] --configfile [config1,config2,config3] --cores N
```

You can still run individual snakemake rules as before. Each of the three steps (explained more below) has its own Snakefile, so you must specify the correct snakefile when running an individual rule.

### Providing externally-generated data

Unlike Tourmaline 1, you can start any of the three workflow steps with data from an external program, so long as it is formatted correctly. For example, if you already have ASV sequences and just want to assign taxonomy with Tourmaline, you can format them for QIIME 2 (code to help with this below) and just provide the file path in your config file.

## Overview

Tourmaline 2 is a modular Snakemake pipeline for processing DNA metabarcoding data. The pipeline consists of three main steps, plus an optional fourth step:

### Step 1. Sequence quality assurance and quality control

* Called "qaqc" in Tourmaline 2 code.
* Processes raw fastq files (paired-end or single-end data).
* Provides sequence quality plots for demultiplexed raw and/or trimmed reads.
* Optionally trims primer sequences from raw reads.
* Creates a QIIME 2 sequence artifact.

See [QA/QC Step](steps/qaqc.md) for details.

### Step 2. Representative sequences (denoising and ASV generation)

* Called "repseqs" in Tourmaline 2 code.
* Generates ASVs using the specified method (DADA2 or Deblur).
* Optional filtering based on length, abundance, and prevalence.
* Produces feature table and representative sequences.

See [Repseqs Step](steps/repseqs.md) for details.

### Step 3. Taxonomy assignment

* Called "taxonomy" in Tourmaline 2 code.
* Generates taxonomic assignments and visualizations.
* Assigns taxonomy using one of four methods:
  * [Naive Bayes classifier as implemented in QIIME 2](https://docs.qiime2.org/2024.10/plugins/available/feature-classifier/classify-sklearn/)
  * [Consensus BLAST as implemented in QIIME 2](https://docs.qiime2.org/2024.10/plugins/available/feature-classifier/classify-consensus-blast/)
  * [Consensus VSEARCH as implemented in QIIME 2](https://docs.qiime2.org/2024.10/plugins/available/feature-classifier/classify-consensus-vsearch/)
  * [Anacapa's Bowtie 2 and BLCA method](https://github.com/limey-bean/Anacapa?tab=readme-ov-file#step-3-taxonomic-assignment-using-bowtie-2-and-blca)

See [Taxonomy Step](steps/taxonomy.md) for details.

### Step 4. Generate bioinformatics metadata

* Creates a file with metadata about the analysis using FAIR eDNA terms.
* File can be read into the [NOAA Ocean DNA Explorer](https://www.ngi.msstate.edu/node).

See [Metadata](metadata.md) for details.

## Quick Start

Tourmaline 2 provides a modular workflow for processing amplicon sequencing data. The pipeline consists of three main steps that can be run together or independently:

1. **QA/QC** - Process raw FASTQ files, optional primer trimming, generate QIIME 2 artifact
2. **Repseqs** - Generate ASVs using DADA2 or Deblur, optional filtering, produce feature table and representative sequences
3. **Taxonomy** - Assign taxonomy using one of four methods, generate visualizations

### Getting Started

1. **Install and Setup**: See [Install and Setup](install.md) for requirements and environment setup.
   - Install QIIME 2 (2024.10) amplicon workflow
   - Create Snakemake conda environment
   - Clone Tourmaline repository

2. **Configuration**: See [Configuration](configuration.md) for setting up config files for each step.
   - Create config files for each step you plan to run
   - Configure parameters for your data type and analysis needs

3. **Running**: See [Running](running.md) for how to use `tourmaline.sh` and examples.
   - Run single steps or all steps together
   - Use the `tourmaline.sh` script for easy execution

### Example Workflow

```bash
# Activate environment
conda activate snakemake-tour2

# Run all steps
./tourmaline.sh -s qaqc,repseqs,taxonomy \
  -c config_01_qaqc.yaml,config_02_repseqs.yaml,config_03_taxonomy.yaml \
  -n 6

# Or run a single step
./tourmaline.sh -s taxonomy -c config_03_taxonomy.yaml -n 6
```

For detailed instructions, see the [Install and Setup](install.md), [Configuration](configuration.md), and [Running](running.md) pages.

## Documentation Structure

- **[Install and Setup](install.md)**: Requirements, conda environments, and getting Tourmaline
- **[Configuration](configuration.md)**: Config file parameters for all three steps
- **[Running](running.md)**: Using `tourmaline.sh` script and examples
- **[Steps](steps/qaqc.md)**: Detailed documentation for each pipeline step
  - [QA/QC](steps/qaqc.md): Sequence quality control and trimming
  - [Repseqs](steps/repseqs.md): ASV generation with DADA2 or Deblur
  - [Taxonomy](steps/taxonomy.md): Taxonomic assignment methods
- **[External Data](external_data.md)**: Providing externally-generated inputs and conversions
- **[Metadata](metadata.md)**: Generating bioinformatics metadata
- **[Troubleshooting](troubleshooting.md)**: Common issues and tips
- **[Citation & Legacy](citation_legacy.md)**: How to cite Tourmaline and v1 resources

## Directory structure

The pipeline creates the following directory structure for outputs:

```
output_dir/
├── [run_name]-qaqc/    # QA/QC outputs (was "samples" in some docs)
├── [run_name]-repseqs/    # Representative sequences outputs
└── [run_name]-taxonomy/   # Taxonomy assignment outputs
```

Each directory contains the relevant outputs for that step of the pipeline.

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
