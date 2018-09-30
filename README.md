<img src="https://upload.wikimedia.org/wikipedia/commons/0/00/Tourmaline-121240.jpg" height=200> <img src="http://melissabessmonroe.com/wp-content/uploads/2014/03/20140303_TourmalineSurfPark128.jpg" height=200>

<!--[![Build Status](https://travis-ci.org/cuttlefishh/tourmaline.svg?branch=master)](https://travis-ci.org/cuttlefishh/tourmaline)

See https://docs.travis-ci.com/user/getting-started/ and https://github.com/biocore/oecophylla/blob/master/.travis.yml for setting up Travis.-->

# tourmaline

Amplicon sequencing is a metagenetics method whereby a single DNA locus in a community of organisms is PCR-amplified and sequenced. Tourmaline is an amplicon sequence processing workflow for Illumina sequence data that uses [QIIME 2](https://qiime2.org) and the software packages it wraps. Tourmaline manages commands, inputs, and outputs using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system.

Two methods of amplicon sequence processing are supported, both of which generate ASVs (amplicon sequence variants, which approximate the "true" or "exact" sequences in a sample) rather than OTUs (operational taxonomic units, which blur sequencing errors and microdiversity through clustering):

* [Deblur](https://github.com/biocore/deblur) is a greedy deconvolution algorithm based on known Illumina read error profiles ([Amir et al., 2017](https://doi.org/10.1128/mSystems.00191-16)).
* [DADA2](https://github.com/benjjneb/dada2) implements a quality-aware model of Illumina amplicon errors to infer sample composition by dividing amplicon reads into partitions consistent with the error model ([Callahan et al., 2016](https://doi.org/10.1038/nmeth.3869)).

Tourmaline is an alternative amplicon pipeline to [Banzai](https://github.com/jimmyodonnell/banzai), which was developed for [MBON](https://github.com/marinebon/MBON) and uses [Swarm](https://github.com/torognes/swarm) for OTU picking.

## Installation

Tourmaline has the following dependencies (installation instructions below):

* Conda
* QIIME 2 version `2018.6` (should work with later versions but has not been tested)
* Snakemake

### QIIME2 and Snakemake

First, if you don't have Conda installed on your machine, install [Miniconda](https://conda.io/miniconda.html) for your operating system (Python 3.7+ version).

Second, install QIIME 2 in a Conda environment, if you haven't already. See the instructions at [qiime2.org](https://docs.qiime2.org/2018.6/install/native/). For example, on macOS these commands will install QIIME 2 inside a Conda environment called `qiime2-2018.6`:

```
wget https://data.qiime2.org/distro/core/qiime2-2018.6-py35-osx-conda.yml
conda env create -n qiime2-2018.6 --file qiime2-2018.6-py35-osx-conda.yml
```

Third, activate your QIIME 2 environment and install Snakemake:

```
source activate qiime2-2018.6
conda install snakemake
```

Finally, clone the Tourmaline repository and rename it to the working directory for your project (replace the directories in ALL CAPS with your directories):

```
cd /PATH/TO
git clone https://github.com/cuttlefishh/tourmaline.git
mv tourmaline PROJECT
```

### Test Data

The Tourmaline repository comes ready to go with test 18S rRNA fastq sequence data and a corresponding reference database.

First you must edit the manifest files `00-data/manifest_se.csv` and `00-data/manifest_pe.csv` to point to the absolute filepaths of the sequences in your local copy of `tourmaline` (which you renamed to `PROJECT`). For example, if the filepath of your project is `/PATH/TO/PROJECT`, these commands will fix the manifest files:

```
cat manifest_pe.csv | sed 's|/Users/luke.thompson/git/tourmaline|/PATH/TO/PROJECT|' > temp
mv temp manifest_pe.csv 
cat manifest_se.csv | sed 's|/Users/luke.thompson/git/tourmaline|/PATH/TO/PROJECT|' > temp
mv temp manifest_se.csv
```

You will also need to edit the configuration file `config.yaml` to reduce the Deblur trim length because the test sequences are only 120 bp, and decrease the subsampling values because the sequencing depth of the test dataset is very low:

```
deblur_trim_length: 100
alpha_max_depth: 50
core_sampling_depth: 50
```

Hint: Before you change `config.yaml`, make a copy called `config_default.yaml` that will stay unchanged. You can always run `diff config_default.yaml config.yaml` to see which parameters you have changed from the defaults.

<!--
THIS DOES NOT FIX DEBLUR WITH TEST DATA, STILL GET IndexError:

deblur_min_reads: 1
deblur_min_size: 1
-->

Now you are ready to test Snakemake. You might start with the DADA2 paired-end workflow:

```
snakemake dada2_pe_denoise
```

You can try any of the DADA2 rules below. *Warning: The current test dataset does not work with Deblur.* Once you are ready to start analyzing your own data, make sure to delete the output directories (`01-imported`, `02-denoised`, `03-repseqs`, `04-diversity`, `05-reports`) generated in the testing process.

### Helper Scripts

Tourmaline comes with the following helper scripts in the `scripts` directory:

* `match_manifest_to_metadata.py` -- Given a metadata/mapping file and a fastq manifest file, generate paired-end and single-end manifest fastq files corresponding to the samples in the metadata/mapping file. This is useful if you want to use Tourmaline, which requires paired-end and/or single-end manifest files, but your data are not demultiplexed (e.g., they are in the "Earth Microbiome Project (EMP) protocol" format for paired-end reads, as described below).
* `detect_amplicon_locus.py` -- Try to guess the amplicon locus of a fasta file based on the first four nucleotides. Used by rule `repseq_detect_amplicon_type`.
* `fastaLengths.pl` -- Calculate the length of every sequence in a fasta file. Used by rule `repseq_length_distribution`.

## Execution

This section describes how to prepare your data and then run the workflow.

### Prepare data

Tourmaline steps covered in this section (logic described below):

* Step 0: Assess and format data

Note: Currently, Tourmaline does not do demultiplexing and quality filtering (support for this may be added later). You should demultiplex and quality filter your fastq data yourself, then gzip the resulting per-sample fastq files.

#### Format metadata

Your metadata, also known as a mapping file, should be a tab-delimited text file (e.g., exported from Excel) with samples as rows and metadata categories as columns. 

Sample metadata should include basic sample information like `collection_timestamp`, `latitude`, and `longitude`, plus categories that describe treatment groups and environmental metadata relevant to your samples. It's important to have rich and complete sample metadata before you begin your analyses.

Processing information should ideally include all of the following columns: `project_name`, `experiment_design_description`, `target_gene`, `target_subfragment`, `pcr_primers`, `pcr_primer_names`, `platform`, `instrument_model`, `run_center`, `run_date`. These will be included in your QC report.

The above columns follow the standards set by [Qiita](https://qiita.ucsd.edu/static/doc/html/gettingstartedguide/index.html). For additional help see [Metadata in QIIME 2](https://docs.qiime2.org/2018.6/tutorials/metadata/), the [EMP Metadata Guide](http://www.earthmicrobiome.org/protocols-and-standards/metadata-guide/), and [QIIMP](https://qiita.ucsd.edu/iframe/?iframe=qiimp) for help formatting your metadata.

#### Format sequence data

Tourmaline supports amplicon sequence data that is already demultiplexed (fastq manifest format). Using the sample names in your metadata file and the absolute filepaths to the forward and reverse demultiplexed sequence files (`.fastq.gz`) for each sample, create a fastq manifest file. See [Fastq Manifest Formats](https://docs.qiime2.org/2018.6/tutorials/importing/#fastq-manifest-formats) (QIIME 2) for instructions for creating this file. If your sequences are not already demultiplexed (e.g., they need to be imported as type `EMPPairedEndSequences`), you can use the commands `qiime tools import` and `qiime demux emp-paired` to demuliplex them, then unzip the archives and merge the manifest files, taking care to change the second column to `absolute-filepath` and ensure the sample IDs match those in your metadata file (the included script `match_manifest_to_metadata.py` can help with this).

Note: While `qiime tools import` supports both `.fastq` and `.fastq.gz` formats, using `.fastq.gz` format is strongly recommended because it is ~5x faster and minimizes disk usage. (Hint: Gzipped files can still be viewed using `zcat` with `less` or `head`.)

#### Set up data directory

We'll use variables to reference our project and databases directories:

```
PRJ=/PATH/TO/PROJECT
DB=/PATH/TO/DATABASES
```

Create a directory for data inside your working directory:

```
cd $PRJ
mkdir 00-data
```

Move your metadata file and fastq manifest file(s) to `$PRJ/00-data`.

#### Create symbolic links to reference sequence and taxonomy data

Reference sequences and taxonomy are used by QIIME 2 to assign taxonomy to the representative sequences you get from Deblur and DADA2. The sequence (full-length and trimmed to primers) and taxonomy files, and their QIIME 2 archive derivatives, can be used for multiple studies and need only exist in one place on your computer. Therefore, it may be best to store these files in a special directory for databases, then create symbolic links in the respective project subdirectories.

Symbolic links to the QIIME 2-formatted reference sequence fasta and taxonomy files can go here (example locations of files shown):

```
cd $PRJ/00-data
ln -s $DB/16s/16s_refseqs.fna refseqs.fna
ln -s $DB/16s/16s_reftax.tsv reftax.tsv
```

The first time you run Tourmaline with a given amplicon locus, the files below will be created. However, for future projects with that same amplicon locus, you can put these files in your databases directory. You will have to make the directory `01-imported` before creating the symbolic links:

```
mkdir $PRJ/01-imported
cd $PRJ/01-imported
ln -s $DB/16s/16s_refseqs.qza refseqs.qza
ln -s $DB/16s/16s_reftax.qza reftax.qza
ln -s $DB/16s/16s_refseqs_515f_806r.qza refseqs_extracted.qza
ln -s $DB/16s/16s_classifier.qza classifier.qza
```

For Snakemake to work with these symbolic links, you may have to run `snakemake --cleanup-metadata <filenames>` on them first.

#### Edit the configfile

The configuration file or `configfile` is `config.yaml`. It must be edited to contain the paths to your data and the parameters you want to use. `Snakefile` and `config.yaml` should describe all the inputs, parameters, and commands needed to produce the desired output.

### Run Snakemake

Tourmaline steps covered in this section (logic described below):

* Step 1: Import data
* Step 2: Denoising
* Step 3: Representative sequence curation
* Step 4: Core diversity analyses
* Step 5: Quality control report

#### Rules

Snakemake works by executing rules, defined in the `Snakefile`. Rules specify commands and outputs but most critically inputs, which dictate which other rules must be run beforehand to generate those inputs. By defining pseudo-rules at the beginning of the `Snakefile`, we can specify desired endpoints as "inputs" that force execution of the whole workflow or just part of it. When a Snakemake command is run, only those rules that need to be executed to produce the requested inputs will be run. Before proceeding, you should familiarize yourself with Snakemake using the [documentation](https://snakemake.readthedocs.io), and create your own simple Snakemake workflow to understand how it works.

Tourmaline provides Snakemake rules for Deblur (single-end) and DADA2 (single-end and paired-end). For each type of processing, the `denoise` rule imports data and runs denoising (steps 1 and 2), the `diversity` rule does representative sequence curation and core diversity analyses (steps 3 and 4), the `stats` rule runs group significance and other tests (optional), and the `report` rule generates the QC report (step 5). Pausing after step 2 allows you to make changes before proceeding:

* Check the table summaries and representative sequence lengths to determine if Deblur or DADA2 parameters need to be modified. If so, you can rename the output directories and then rerun the `denoise` rule.
* View the table visualization to decide an appropriate subsampling (rarefaction) depth. Then modify the parameters `alpha_max_depth` and `core_sampling_depth` in `config.yaml`.
* Filter your biom table and representative sequences to remove unwanted sequences. For example, if your amplicon is 16S rRNA, you may want to filter out chloroplast/mitochondria sequences. You should keep the same filenames so that Snakemake will recognize them; you can save the old versions with different names if you don't want to overwrite them.

##### Deblur (single-end)

```
# steps 1-2
snakemake deblur_se_denoise

# steps 3-4
snakemake deblur_se_diversity

# statistical analyses
snakemake deblur_se_stats

# step 5
snakemake deblur_se_report
```

##### DADA2 (single-end)

```
# steps 1-2
snakemake dada2_se_denoise

# steps 3-4
snakemake dada2_se_diversity

# statistical analyses
snakemake dada2_se_stats

# step 5
snakemake dada2_se_report
```

##### DADA2 (paired-end)

```
# steps 1-2
snakemake dada2_pe_denoise

# steps 3-4
snakemake dada2_pe_diversity

# statistical analyses
snakemake dada2_pe_stats

# step 5
snakemake dada2_pe_report
```

That's it. Just run one of these commands and let Snakemake do its magic. The results will be placed in organized directories inside your working directory:

```
01-imported
02-denoised
03-repseqs
04-diversity
05-reports
```

## Logic

In plain English, this is the logic behind Tourmaline. Starting with Step 1, these steps correspond to the rules (commands) in the `Snakefile`.

### Step 0: Assess and format data

Answer the following questions to determine the best parameters for processing and to be able to evaluate the success of your completed workflow.

#### Assess amplicon locus

* What is the locus being amplified, and what are the primer sequences?
* How much sequence variation is expected for this locus (and primer sites) and dataset?
* Is the expected sequence variation enough to answer my question?
* What is the expected amplicon size for this locus and dataset?

#### Assess sequence data

* What type and length of sequencing was used? (e.g., MiSeq 2x150bp)
* Do I have long enough sequening to do paired-end analysis, or do I have to do single-end analysis only?
* What sequence pre-processing has been done already: Demultiplexing? Quality filtering and FastQC? Primer removal? Merging of paired reads?

#### Assess sample set and metadata

* Is my metadata file complete? Are the relevant parameters of my dataset present as numeric or categorical variables?
* Do I have enough samples in each group of key metadata categories to determine an effect?

#### Format metadata and sequence data

* Is my metadata file properly formatted? See [Metadata in QIIME 2](https://docs.qiime2.org/2018.6/tutorials/metadata/), the [EMP Metadata Guide](http://www.earthmicrobiome.org/protocols-and-standards/metadata-guide/), and [QIIMP](https://qiita.ucsd.edu/iframe/?iframe=qiimp) for help formatting your metadata.
* Is my sequence data demultiplexed, in `.fastq.gz` format, and described in a QIIME 2 fastq manifest file? See [Fastq Manifest Formats](https://docs.qiime2.org/2018.6/tutorials/importing/#fastq-manifest-formats) from QIIME 2 for instructions for creating this file.
* Are my reference sequences and taxonomy properly formatted for QIIME 2?
* Is my config file updated with the file paths and parameters I want to use?

### Step 1: Import data

Import sequence data and metadata for QIIME 2 processing. It is assumed that any preprocessing and formatting of data have already been done.

* Amplicon sequence data: import into QIIME 2 artifact.
* Reference sequence data: import into QIIME 2 artifact.
* Metadata: import into QIIME 2 artifact.

### Step 2: Denoising

Run Deblur or DADA2 to generate ASV feature tables and representative sequences. This is akin to OTU picking.

* Denoise amplicon data using Deblur or DADA2 to generate ASV feature tables (BIOM).
* Use paired-end mode (DADA2 only) if sequence and amplicon lengths permit.

Perform quality control of output tables and representative sequences. Manually evaluate output files to determine subsampling (rarefaction) depth, decide if denoising parameters need to be changed, or if output need to be filtered. QC results will also go into the QC report later.

* Summarize tables and representative sequences (ASVs).
* Manually view table summary to determine appropriate rarefaction depth.
* Manually view representative sequences summary to determine length distribution of sequences.

### Step 3: Representative sequence curation

Generate a phylogenetic tree of ASV sequences, and identify the taxonomy (phylum, class, order, etc.) of each ASV.

* Build a phylogenetic tree of ASV sequences, or insert ASV sequences into an existing tree for your amplicon locus.
* Assign taxonomy to ASVs using a reference database for your amplicon locus.

### Step 4: Core diversity analyses

First consult table summary and run alpha rarefaction to decide on a rarefaction depth. Then do the major alpha/beta diversity analyses and taxonomy summary.

* Alpha diversity: alpha rarefaction, diversity metrics (evenness, Shannon, Faith's PD, observed sequences), alpha group significance.
* Beta diversity: distance matrices (un/weighted UniFrac, Bray-Curtis, Jaccard), principal coordinates, Emperor plots, beta group significance.
* Taxonomy barplots.

### Step 4.1: Statistical analyses

Run statistical tests on your data, such as group significance tests. These rules are optional, and the output does not go into the QC report. (This functionality is currently minimal.)

### Step 5: Quality control report

After completing processing and core analyses, determine if the results make sense.

* How many sequences did I start with, and how many are left after denoising?
* Are the representative sequences of similar length or of very different lengths?
* Do the sequence alignment and tree look reasonable?
* Do samples cluster in an expected way in PCoA space?
* Do the taxonomic profiles match expected taxonomic compositions?
