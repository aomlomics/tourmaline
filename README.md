<img src="tourmaline_banner.png" alt="tourmaline_banner" width="100%"/>

## Tourmaline

Tourmaline is an amplicon (metabarcoding) sequence processing workflow for Illumina sequence data. It uses the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system as a wrapper for [QIIME 2](https://qiime2.org) and additional shell and Python scripts.

### Why should I use Tourmaline?

* **QIIME 2.** The core commands of Tourmaline, including the [DADA2](https://benjjneb.github.io/dada2/index.html) package, are all commands of QIIME 2, one of the most popular amplicon sequence analysis software tools available. You can print all of the QIIME 2 and other shell commands of your workflow before or while running the workflow.
* **Snakemake.** Managing the workflow with Snakemake provides several benefits: 
  - **Configuration file** contains all parameters in one file, so you can see what your workflow is doing and make changes for a subsequent run.
  - **Directory structure** is the same for every Tourmaline run, so you always know where your outputs are.
  - **On-demand commands** mean that only the commands required for output files not yet generated are run, saving time and computation when re-running part of a workflow.
* **Parameter optimization.** The configuration file and standard directory structure make it simple to test and compare different parameter sets to optimize your workflow. Included code helps choose read truncation parameters and identify outliers in representative sequences (ASVs).
* **Reports.** Every Tourmaline run produces an HTML report containing a summary of your metadata and outputs, with links to web-viewable QIIME 2 visualization files.
* **Tourmaline Toolkit.** Analyze multiple outputs programmatically using the provided code and notebooks written in R and Python.

Ready to get started? Visit the [Wiki](https://github.com/lukenoaa/tourmaline/wiki) for a detailed guide on using Tourmaline. If you're feeling bold and want to get started right away, check out the Quick Start instructions below.

## Quick Start



## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

