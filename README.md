<img src="png/tourmaline_banner.png" alt="png/tourmaline_banner" width="100%"/>

<img src="png/figure1.png" alt="png/figure1" width="100%"/>

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
* **Downstream analysis.** Analyze the output of single or multiple Tourmaline runs programmatically, with qiime2R in R or the QIIME 2 Artifact API in Python, using the provided R and Python notebooks or your own code.

Ready to get started? Visit the [Wiki](https://github.com/lukenoaa/tourmaline/wiki) for a detailed guide on using Tourmaline. If you're feeling bold and want to get started right away, check out the Quick Start instructions below.

## Quick Start

Tourmaline provides Snakemake rules for DADA2 (single-end and paired-end) and Deblur (single-end). For each type of processing, there are four steps:

1. the *denoise* rule imports FASTQ data and runs denoising, generating a feature table and representative sequences;
2. the *taxonomy* rule assigns taxonomy to representative sequences;
3. the *diversity* rule does representative sequence curation, core diversity analyses, and alpha and beta group significance; and
4. the *report* rule generates an HTML report of the outputs plus metadata, inputs, and parameters.

Steps 2–4 have *unfiltered* and *filtered* modes, the difference being that in the *taxonomy* step of *filtered* mode, undesired taxonomic groups or individual sequences from the representative sequences and feature table are removed. The *diversity* and *report* rules are the same for *unfiltered* and *filtered* modes, except the output goes into separate subdirectories.

### Install

Tourmaline requires a Conda installation of QIIME 2, Snakemake, and other dependencies:

```
wget https://data.qiime2.org/distro/core/qiime2-2020.8-py36-osx-conda.yml
conda env create -n qiime2-2020.8 --file qiime2-2020.8-py36-osx-conda.yml
conda activate qiime2-2020.8
conda install -c bioconda snakemake biopython tabulate pandoc tabview
conda install -c bioconda bioconductor-msa bioconductor-odseq
pip install git+https://github.com/biocore/empress.git
qiime dev refresh-cache
```

Alternatively, use the Docker image to create a container with all the dependencies:

```
docker pull aomlomics/tourmaline
docker run -it aomlomics/tourmaline
```

### Setup

If this is your first time running Tourmaline, you'll need to set up your directory. Simplified instructions are below, but see the Wiki's [Setup](https://github.com/lukenoaa/tourmaline/wiki/3-Setup) page for complete instructions. 

Start by cloning the Tourmaline directory and files:

```
git clone https://github.com/aomlomics/tourmaline.git
```

#### Setup for the test data

The test data (16 samples of paired-end 16S rRNA data with 1000 sequences per sample) comes with your cloned copy of Tourmaline. It's fast to run and will verify that you can run the workflow.

Download reference database sequence and taxonomy files, named `refseqs.qza` and `reftax.qza` (QIIME 2 archives), in `01-imported`:

```
cd tourmaline/01-imported
wget https://data.qiime2.org/2020.8/common/silva-138-99-seqs-515-806.qza
wget https://data.qiime2.org/2020.8/common/silva-138-99-tax-515-806.qza
ln -s silva-138-99-seqs-515-806.qza refseqs.qza
ln -s silva-138-99-tax-515-806.qza reftax.qza
```

Edit FASTQ manifests `manifest_se.csv` and `manifest_pe.csv` in `00-data` so file paths match the location of your `tourmaline` directory. In the command below, replace `/path/to` with the location of your `tourmaline` directory—or skip this step if you are using the Docker container and you cloned `tourmaline` into `/data`:

```
cd ../00-data
cat manifest_pe.csv | sed 's|/data/tourmaline|/path/to/tourmaline|' > temp
mv temp manifest_pe.csv 
cat manifest_pe.csv | grep -v "reverse" > manifest_se.csv
```

Create a symbolic link from `Snakefile_mac` or `Snakefile_linux` to `Snakefile`, depending on your operating system (if you are using the Docker container, your operating system is Linux):

```
cd ..
ln -s Snakefile_linux Snakefile
```

Go to **Run Snakemake**.

#### Setup for your data

If you're ready to run your own data, the setup is similar to what you did for the test data:

* Put reference database taxonomy and FASTA files in `00-data` or imported QIIME 2 archives in `01-imported`.
* Edit FASTQ manifests `manifest_se.csv` and `manifest_pe.csv` so file paths point to your .fastq.gz files (they can be anywhere on your computer) and sample names match the metadata file.
* Edit metadata file `metadata.tsv` to contain your sample names and any relevant metadata for your samples.
* Edit configuration file `config.yaml` to change PCR locus/primers, DADA2/Deblur parameters, and rarefaction depth.
* Create a symbolic link from `Snakefile_mac` or `Snakefile_linux` (depending on your system) to `Snakefile`.
* Go to **Run Snakemake**.

Note: If you've run Tourmaline on your dataset before, you can skip the steps above and initialize a new Tourmaline directory (e.g., `tourmaline-new`) with the files and symlinks of the existing one (e.g., `tourmaline-existing`) using the command below:

```
cd /path/to/tourmaline-new
scripts/initialize_dir_from_existing_tourmaline_dir.sh /path/to/tourmaline-existing
```

Just remember to make any changes to your configuration file before you run Snakemake.

### Run Snakemake

Shown here is the DADA2 paired-end workflow. See the Wiki's [Run](https://github.com/lukenoaa/tourmaline/wiki/4-Run) page for complete instructions on all steps, denoising methods, and filtering modes.

From the `tourmaline` directory (which you may rename), run Snakemake with the *denoise* rule as the target:

```
snakemake dada2_pe_denoise
```

Pausing after the *denoise* step allows you to make changes before proceeding:

* Check the table summaries and representative sequence lengths to determine if DADA2 or Deblur parameters need to be modified. If so, you can rename or delete the output directories and then rerun the *denoise* rule.
* View the table visualization to decide an appropriate subsampling (rarefaction) depth. Then modify the parameters "alpha_max_depth" and "core_sampling_depth" in `config.yaml`.
* Decide whether to filter your feature table and representative sequences by taxonomy or feature ID. After the *taxonomy* step, you can examine the taxonomy summary and bar plot to aid your decision. If you do filter your data, all output from that point on will go in a separate folder so you can compare output with and without filtering.

#### Unfiltered mode

Continue the workflow without filtering (for now). If you are satisfied with your parameters and files, run the *taxonomy* rule (for unfiltered data):

```
snakemake dada2_pe_taxonomy_unfiltered
```

Next, run the *diversity* rule (for unfiltered data):

```
snakemake dada2_pe_diversity_unfiltered
```

Finally, run the *report* rule (for unfiltered data):

```
snakemake dada2_pe_report_unfiltered
```

#### Filtered mode

Filtering is done on representative sequences and the feature table, and downstream outputs will be filtered; the taxonomy file itself is not filtered. Filtering can be done by taxonomy keywords and/or by feature IDs of specific representative sequences. To filter by these two methods, before running *taxonomy_filtered*, do this:

* Taxonomy keyword: Place the keywords in `config.yaml` in the field "exclude_terms", separated by commas. Searching is not case-sensitive.
* Specific representative sequences: Go to `2-output-dada2-pe-unfiltered/02-alignment-tree` and copy or merge `repseqs_to_filter_outliers.tsv` and/or `repseqs_to_filter_unassigned.tsv` to  `00-data/repseqs_to_filter_dada2-pe.tsv`. If merging the two files, take care to remove duplicate feature IDs, because duplicates will cause the filtering step to fail.

Now we are ready to filter the representative sequences and feature table, generate new summaries, and generate a new taxonomy bar plot, by running the *taxonomy* rule (for filtered data):

```
snakemake dada2_pe_taxonomy_filtered
```

Next, run the *diversity* rule (for filtered data):

```
snakemake dada2_pe_diversity_filtered
```

Finally, run the *report* rule (for filtered data):

```
snakemake dada2_pe_report_filtered
```

#### Troubleshooting

* The whole workflow should take ~3–5 minutes to complete with the test data. A normal dataset may take several hours to complete.
* If any of the above commands don't work, read the error messages carefully, try to figure out what went wrong, and attempt to fix the offending file. A common issue is the file paths in your FASTQ manifest file need to be updated.
* Do not use the `--cores` option. Tourmaline should be run with 1 core (default).
* If you are running in a Docker container and you get an error like "Signals.SIGKILL: 9", you probably need to give Docker more memory. See the Wiki section on [Installation options](https://github.com/lukenoaa/tourmaline/wiki/2-Install#installation-options).

#### Power tips

* The whole workflow can be run with just the command `snakemake dada2_pe_report_unfiltered`  (without filtering representative sequences) or  `snakemake dada2_pe_report_filtered`  (after filtering representative sequences). Warning: If your parameters are not optimized, the results will be suboptimal (garbage in, garbage out).
* If you want to make a fresh run and not save the previous output, simply delete the output directories (e.g., `02-output-{method}-{filter}` and `03-report`) generated in the previous run.

## License

Software code created by U.S. Government employees is not subject to copyright in the United States (17 U.S.C. §105). The United States/Department of Commerce reserve all rights to seek and obtain copyright protection in countries other than the United States for Software authored in its entirety by the Department of Commerce. To this end, the Department of Commerce hereby grants to Recipient a royalty-free, nonexclusive license to use, copy, and create derivative works of the Software outside of the United States.

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
