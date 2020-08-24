<img src="tourmaline_banner.png" alt="tourmaline_banner" width="100%"/>

## Tourmaline

Tourmaline is an amplicon (metabarcoding) sequence processing workflow for Illumina sequence data. It uses the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system as a wrapper for [QIIME 2](https://qiime2.org) and additional shell and Python scripts.

### Why should I use Tourmaline?

* **QIIME 2.** The core commands of Tourmaline are all commands of QIIME 2, one of the most popular amplicon sequence analysis software tools availalbe. You can print all of the QIIME 2 and other shell commands of your workflow before or while running the workflow.
* **Snakemake.** Managing the workflow with Snakemake provides several benefits: 
  - **Configuration file** contains all parameters in one file, so you can see what your workflow is doing and make changes for a subsequent run.
  - **Directory structure** is the same for every Tourmaline run, so you always know where your outputs are.
  - **On-demand commands** mean that only the commands required for output files not yet generated are run, saving time and computation when re-running part of a workflow.
* **Parameter optimization.** The configuaration file and standard directory structure make it simple to test and compare different parameter sets to optimize your workflow.
* **Reports.** Every Tourmaline run produces an HTML report containing a summary of your metadata and outputs, with links to web-viewable QIIME 2 visualization files.
* **Tourmaline Toolkit.** Analyze multiple outputs programmatically using the provided code and notebooks written in R and Python.

Ready to get started? Visit the [Wiki](https://github.com/lukenoaa/tourmaline/wiki) for a detailed guide on using Tourmaline. If you're feeling bold and want to get started right away, check out the Quick Start instructions below.

## Quick Start

Tourmaline provides Snakemake rules for DADA2 (single-end and paired-end) and Deblur (single-end). For each type of processing, the "denoise" rule imports data and runs denoising; the "diversity" rule does representative sequence curation, core diversity analyses, and alpha and beta group significance; and the "report" rule generates the QC report. 

### Setup

Start by cloning the Tourmaline directory and files:

```
git clone https://github.com/aomlomics/tourmaline.git
```

If this is your first time running Tourmaline, you'll need to set up your directory. See the Wiki's [Setup](https://github.com/lukenoaa/tourmaline/wiki/3-Setup) page for instructions. Briefly, to process the **Test data**:

* Put reference database taxonomy and FASTA (as imported QIIME 2 archives) in `01-imported`.
* Edit FASTQ manifests `manifest_se.csv` and `manifest_pe.csv` in `00-data` so file paths match the location of your `tourmaline` directory.
* Create a symbolic link from `Snakefile_mac` or `Snakefile_linux` (depending on your system) to `Snakefile`.

Or to process **Your data**:

* Put reference database taxonomy and FASTA files in `00-data` or imported QIIME 2 archives in `01-imported`.
* Edit FASTQ manifests `manifest_se.csv` and `manifest_pe.csv` so file paths point to your .fastq.gz files (they can be anywhere on your computer) and sample names match the metadata file.
* Edit metadata file `metadata.tsv` to contain your sample names and any relevant metadata for your samples.
* Edit configuration file `config.yaml` to change PCR locus/primers, DADA2/Deblur parameters, and rarefaction depth.
* Create a symbolic link from `Snakefile_mac` or `Snakefile_linux` (depending on your system) to `Snakefile`.

If you've run Tourmaline on your dataset before, you can initialize a new Tourmaline directory with the files and symlinks of an existing one using the command below:

```
cd /PATH/TO/NEW/TOURMALINE
scripts/initialize_dir_from_existing_tourmaline_dir.sh /PATH/TO/EXISTING/TOURMALINE
# then make any changes to your configuration before running
```

### Run Snakemake

Shown here is the DADA2 paired-end workflow. From the `tourmaline` directory (which you may rename), run Snakemake with the "denoise" rule as the target:

```
snakemake dada2_pe_denoise
```

Pausing after the "denoise" step allows you to make changes before proceeding:

* Check the table summaries and representative sequence lengths to determine if DADA2 or Deblur parameters need to be modified. If so, you can rename or delete the output directories and then rerun the "denoise" rule.
* View the table visualization to decide an appropriate subsampling (rarefaction) depth. Then modify the parameters "alpha_max_depth" and "core_sampling_depth" in `config.yaml`.
* Filter your biom table and representative sequences to remove unwanted sequences. For example, if your amplicon is 16S rRNA, you may want to filter out chloroplast/mitochondria sequences. You should keep the same filenames so that Snakemake will recognize them; you can save the old versions with different names if you don't want to overwrite them.

After you are satisfied with your parameters and files, run the "diversity" rule:

```
snakemake dada2_pe_diversity
```

Finally, run the "report" rule:

```
snakemake dada2_pe_report
```

If any of the above commands don't work, read the error messages carefully, try to figure out what went wrong, and attempt to fix the offending file. A common issue is the file paths in your FASTQ manifest file need to be updated.

If you want to make a fresh run and not save the previous output, simply delete the output directories (e.g., `02-output-{method}` and `03-report`) generated in the previous run.

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

