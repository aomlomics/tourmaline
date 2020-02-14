<img src="https://upload.wikimedia.org/wikipedia/commons/0/00/Tourmaline-121240.jpg" height=200> <img src="http://melissabessmonroe.com/wp-content/uploads/2014/03/20140303_TourmalineSurfPark128.jpg" height=200>

<!--[![Build Status](https://travis-ci.org/cuttlefishh/tourmaline.svg?branch=master)](https://travis-ci.org/cuttlefishh/tourmaline)

See https://docs.travis-ci.com/user/getting-started/ and https://github.com/biocore/oecophylla/blob/master/.travis.yml for setting up Travis.-->

# tourmaline

Amplicon sequencing is a metagenetics method whereby a single DNA locus in a community of organisms is PCR-amplified and sequenced. Tourmaline is an amplicon sequence processing workflow for Illumina sequence data that uses [QIIME 2](https://qiime2.org) and the software packages it wraps. Tourmaline manages commands, inputs, and outputs using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system.

#### Amplicon sequence variants

Two methods of amplicon sequence processing are supported, both of which generate ASVs (amplicon sequence variants, which approximate the "true" or "exact" sequences in a sample) rather than OTUs (operational taxonomic units, which blur sequencing errors and microdiversity through clustering):

* [Deblur](https://github.com/biocore/deblur) is a greedy deconvolution algorithm based on known Illumina read error profiles ([Amir et al., 2017](https://doi.org/10.1128/mSystems.00191-16)).
* [DADA2](https://github.com/benjjneb/dada2) implements a quality-aware model of Illumina amplicon errors to infer sample composition by dividing amplicon reads into partitions consistent with the error model ([Callahan et al., 2016](https://doi.org/10.1038/nmeth.3869)).

#### QIIME 2 and Snakemake commands

QIIME 2 shell commands are provided for reference in [`commands.txt`](https://github.com/NOAA-AOML/tourmaline/blob/master/commands.txt). Exact commands executed by the workflow are in [`Snakefile`](https://github.com/NOAA-AOML/tourmaline/blob/master/Snakefile). With the help of Snakemake, Tourmaline provides rapid and reproducible workflows for processing amplicon sequence data, storing output files in a logical directory structure.

#### Ready for MBON

Tourmaline is an alternative amplicon 'pipeline' to [Banzai](https://github.com/jimmyodonnell/banzai), which was developed for [MBON](https://github.com/marinebon/MBON) (Marine Biodiversity Observation Network). Banzai uses [Swarm](https://github.com/torognes/swarm) for OTU picking. Tourmaline supports both [Deblur](https://github.com/biocore/deblur) and [DADA2](https://github.com/benjjneb/dada2) for denoising to generate ASVs. In the future, Tourmaline could be extended to other OTU/ASV picking algorithms if they are added to QIIME 2.

## Installation

Tourmaline requires the following software:

* Conda
* QIIME 2 version 2019.7
* Snakemake

### Conda

First, if you don't have Conda installed on your machine, install [Miniconda](https://conda.io/miniconda.html) for your operating system (Python 3.7+ version).

### QIIME 2

Second, install QIIME 2 in a Conda environment, if you haven't already. See the instructions at [qiime2.org](https://docs.qiime2.org/2019.7/install/native/). For example, on macOS these commands will install QIIME 2 inside a Conda environment called `qiime2-2019.7` (for Linux, change "osx" to "linux"):

```
wget https://data.qiime2.org/distro/core/qiime2-2019.7-py36-osx-conda.yml
conda env create -n qiime2-2019.7 --file qiime2-2019.7-py36-osx-conda.yml
```

### Snakemake

Third, activate your QIIME 2 environment and install Snakemake:

```
conda activate qiime2-2019.7
conda install -c bioconda snakemake
```

Finally, you will install Tourmaline. "Installation" here is really just copying the files to your computer. You will do this in the next step by "cloning" the GitHub repository.

## Setup

### Clone the Tourmaline Repository

Navigate to your project directory and clone the Tourmaline repository there. In the example below and following steps, replace "/PATH/TO/PROJECT" with the full path to your project directory (e.g., "$HOME/workshop-2019.11"):

```
cd $HOME/workshop-2019.11
git clone https://github.com/NOAA-AOML/tourmaline.git
```

You might want to rename the directory `tourmaline` to something else before running the test data, for example (hint: you can do this with your projects to have different copies of Tourmaline with different sample sets or databases):

```
mv tourmaline tourmaline-test
```

Now change directories to `tourmaline-test`:

```
cd tourmaline-test
```

### Snakefile: Mac or Linux

Tourmaline comes with support for both Mac and Linux operating systems. Due to some difference in the way similar commands work on the two systems (e.g., md5/md5sum, gzcat/zcat), two different Snakefiles are provided. Activate the Snakefile for your system by setting a symbolic link:

```
# if you have a mac system
ln -s Snakefile_mac Snakefile
# if you have a linux system
ln -s Snakefile_linux Snakefile
```

You can see where your symbolic links point with the command `ls -l`. When you run Snakemake/Tourmaline, it will use whichever file `Snakefile` points to.

### Test Data

The Tourmaline repository comes ready to go with test 18S rRNA fastq sequence data and a corresponding reference database. 

To run the test data, you must edit the manifest files `00-data/manifest_se.csv` and `00-data/manifest_pe.csv` to point to the absolute filepaths of the sequences in your local copy of `tourmaline` (which you renamed to `tourmaline-test`). For example, if the filepath of your project is `$HOME/workshop-2019.11`, these commands will fix the manifest files (change `$HOME` to the absolute path of your home directory):

```
cd $HOME/workshop-2019.11/tourmaline-test/00-data
cat manifest_pe.csv | sed 's|/PATH/TO/PROJECT/tourmaline|$HOME/workshop-2019.11/tourmaline-test|' > temp
mv temp manifest_pe.csv 
cat manifest_se.csv | sed 's|/PATH/TO/PROJECT/tourmaline|$HOME/workshop-2019.11/tourmaline-test|' > temp
mv temp manifest_se.csv
```

You will need to edit the configuration file `config.yaml` to decrease the subsampling values because the sequencing depth of the test dataset is very low, and if you are testing Deblur reduce the Deblur trim length because the test sequences are only 120 bp in length:

```
deblur_trim_length: 100
...
alpha_max_depth: 50
core_sampling_depth: 50
```

Hint: Before you change `config.yaml`, make a copy called `config_default.yaml` that will stay unchanged. You can always run `diff config_default.yaml config.yaml` to see which parameters you have changed from the defaults.

### Run A Test

Now you are ready to test Snakemake. You might start with the DADA2 paired-end workflow. From your directory `$HOME/workshop-2019.11/tourmaline-test`, run Snakemake with the Denoise rule as the target:

```
snakemake dada2_pe_denoise
```

If that works, try the next target rule in the DADA2 paired-end workflow, the Diversity rule:

```
snakemake dada2_pe_diversity
```

Then try the Stats rule:

```
snakemake dada2_pe_stats
```

Finally try the Report rule:

```
snakemake dada2_pe_report
```

If any of the above commands don't work, read the error messages carefully, try to figure out what went wrong, and attempt to fix the offending file (often the file paths in your fastq manifest need to be changed). Example output from the above commands is found in the directory `example-output`.

You can also try DADA2 single-end. Deblur (command `snakemake deblur_se_denoise`) currently produces an error with the test data, but it should work with normal experimental data. 

<!--
BELOW DOES NOT FIX DEBLUR WITH TEST DATA -- STILL PRODUCES ERROR: IndexError:

deblur_min_reads: 1
deblur_min_size: 1
-->

If you want to use `tourmaline-test` to analyze your own data after testing, make sure to delete the output directories (`01-imported`, `02-denoised`, `03-repseqs`, `04-diversity`, `05-reports`) generated in the testing process.

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

The above columns follow the standards set by [Qiita](https://qiita.ucsd.edu/static/doc/html/gettingstartedguide/index.html). For additional help see [Metadata in QIIME 2](https://docs.qiime2.org/2019.7/tutorials/metadata/), the [EMP Metadata Guide](http://www.earthmicrobiome.org/protocols-and-standards/metadata-guide/), and [QIIMP](https://qiita.ucsd.edu/iframe/?iframe=qiimp) for help formatting your metadata.

#### Format sequence data

Tourmaline supports amplicon sequence data that is already demultiplexed (fastq manifest format). Using the sample names in your metadata file and the absolute filepaths to the forward and reverse demultiplexed sequence files (`.fastq.gz`) for each sample, create a fastq manifest file. See [Fastq Manifest Formats](https://docs.qiime2.org/2019.7/tutorials/importing/#fastq-manifest-formats) (QIIME 2) for instructions for creating this file. If your sequences are not already demultiplexed (e.g., they need to be imported as type `EMPPairedEndSequences`), you can use the commands `qiime tools import` and `qiime demux emp-paired` to demuliplex them, then unzip the archives and merge the manifest files, taking care to change the second column to `absolute-filepath` and ensure the sample IDs match those in your metadata file (the included script `match_manifest_to_metadata.py` can help with this).

Note: While `qiime tools import` supports both `.fastq` and `.fastq.gz` formats, using `.fastq.gz` format is strongly recommended because it is $HOME5x faster and minimizes disk usage. (Hint: Gzipped files can still be viewed using `zcat` with `less` or `head`.)

#### Set up data directory

We'll use variables to reference our project and databases directories:

```
PRJ=/PATH/TO/PROJECT
DB=/PATH/TO/DATABASES
```

Clone a new copy of the Tourmaline repository inside your working directory:

```
cd $PRJ
git clone https://github.com/NOAA-AOML/tourmaline.git
cd tourmaline
```

#### Copy metadata and manifest files or create symbolic links to them

Copy your metadata file and fastq manifest file(s) to `$PRJ/tourmaline/00-data` or create symbolic links to them. Use the code below to remove the test files and create links to your own metadata and fastq manifest files (the files in $PRJ/metadata can be named however you like):

```
cd $PRJ/tourmaline/00-data
rm metadata.tsv manifest_se.csv manifest_pe.csv
ln -s $PRJ/metadata/metadata.tsv metadata.tsv
ln -s $PRJ/metadata/manifest_se.csv manifest_se.csv
ln -s $PRJ/metadata/manifest_pe.csv manifest_pe.csv
```

#### Create symbolic links to reference sequence and taxonomy data

Reference sequences and taxonomy are used by QIIME 2 to assign taxonomy to the representative sequences you get from Deblur and DADA2. The sequence (full-length and trimmed to primers) and taxonomy files, and their QIIME 2 archive derivatives, can be used for multiple studies and need only exist in one place on your computer. Therefore, it may be best to store these files in a special directory for databases, then create symbolic links in the respective project subdirectories.

Symbolic links to the QIIME 2-formatted reference sequence fasta and taxonomy files can go here (example locations of files shown):

```
cd $PRJ/tourmaline/00-data
rm refseqs.fna reftax.tsv
ln -s $DB/16s/16s_refseqs.fna refseqs.fna
ln -s $DB/16s/16s_reftax.tsv reftax.tsv
```

The first time you run Tourmaline with a given amplicon locus, the files below will be created. However, for future projects with that same amplicon locus, you can put these files in your databases directory. You will have to make the directory `01-imported` before creating the symbolic links:

```
mkdir $PRJ/tourmaline/01-imported
cd $PRJ/tourmaline/01-imported
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

Snakemake works by executing rules, defined in the `Snakefile`. Rules specify commands and outputs but most critically inputs, which dictate which other rules must be run beforehand to generate those inputs. By defining pseudo-rules at the beginning of the `Snakefile`, we can specify desired endpoint targets as "inputs" that force execution of the whole workflow or just part of it. When a Snakemake command is run, only those rules that need to be executed to produce the requested target will be run. Before proceeding, you should familiarize yourself with Snakemake using the [documentation](https://snakemake.readthedocs.io) or follow the tutorial [here](https://github.com/cuttlefishh/tutorials/tree/master/snakemake) to create your own simple Snakemake workflow and understand how it works.

#### Directed Acyclic Graph (DAG)

Snakemake provides the command option `--dag` to generate a directed acyclic graph (DAG) of the jobs (rules) that will be run. The DAG is basically a graph that shows the flow and order of rules to reach your desired target. For example, to get a graph (PNG file) of the rules run when setting the rule `snakemake dada2_pe_denoise` as your target, run this command:

```
snakemake dada2_pe_denoise --dag | dot -Tpng > dag.png
```

#### Tourmaline Rules

Tourmaline provides Snakemake rules for Deblur (single-end) and DADA2 (single-end and paired-end). For each type of processing, the `denoise` rule imports data and runs denoising (steps 1 and 2), the `diversity` rule does representative sequence curation and core diversity analyses (steps 3 and 4), the `stats` rule runs group significance and other tests (optional), and the `report` rule generates the QC report (step 5). Pausing after step 2 allows you to make changes before proceeding:

* Check the table summaries and representative sequence lengths to determine if Deblur or DADA2 parameters need to be modified. If so, you can rename the output directories and then rerun the `denoise` rule.
* View the table visualization to decide an appropriate subsampling (rarefaction) depth. Then modify the parameters `alpha_max_depth` and `core_sampling_depth` in `config.yaml`.
* Filter your biom table and representative sequences to remove unwanted sequences. For example, if your amplicon is 16S rRNA, you may want to filter out chloroplast/mitochondria sequences. You should keep the same filenames so that Snakemake will recognize them; you can save the old versions with different names if you don't want to overwrite them.

##### DADA2 (paired-end)

```
# steps 1-2
snakemake dada2_pe_denoise

# steps 3-4
snakemake dada2_pe_diversity

# step 4.1
snakemake dada2_pe_stats

# step 5
snakemake dada2_pe_report
```

##### DADA2 (single-end)

```
# steps 1-2
snakemake dada2_se_denoise

# steps 3-4
snakemake dada2_se_diversity

# step 4.1
snakemake dada2_se_stats

# step 5
snakemake dada2_se_report
```

##### Deblur (single-end)

```
# steps 1-2
snakemake deblur_se_denoise

# steps 3-4
snakemake deblur_se_diversity

# step 4.1
snakemake deblur_se_stats

# step 5
snakemake deblur_se_report
```

That's it. Just run one of these commands and let Snakemake do its magic.

## Output

The results will be placed in organized directories inside your working directory:

```
01-imported
02-denoised
03-repseqs
04-diversity
05-reports
```

The output files of each command (shown for DADA2 paired-end) are as follows:

##### dada2_pe_denoise (steps 1-2)

```
01-imported/fastq_pe.qza
01-imported/fastq_illumina_run.txt
01-imported/fastq_illumina_run.log
02-denoised/data2-pe/stats.qza
02-denoised/data2-pe/table.qza
02-denoised/data2-pe/representative_sequences.qza
02-denoised/data2-pe/table.qzv
02-denoised/data2-pe/table.biom
02-denoised/data2-pe/table_summary_features.txt
02-denoised/data2-pe/table_summary_samples.txt
02-denoised/data2-pe/representative_sequences.fasta
02-denoised/data2-pe/representative_sequences.qzv
02-denoised/data2-pe/representative_sequences_amplicon_type.txt
02-denoised/data2-pe/representative_sequences_lengths.txt
02-denoised/data2-pe/representative_sequences_lengths_describe.tsv
02-denoised/dada2-pe/representative_sequences_md5_status.txt
```

##### dada2_pe_diversity (steps 3-4)

```
01-imported/refseqs.qza
01-imported/reftax.qza
01-imported/classifier.qza
03-repseqs/dada2-pe/taxonomy.qza
03-repseqs/dada2-pe/taxonomy.qzv
03-repseqs/dada2-pe/aligned_representative_sequences.qza
03-repseqs/dada2-pe/masked_aligned_representative_sequences.qza
03-repseqs/dada2-pe/unrooted_tree.qza
03-repseqs/dada2-pe/rooted_tree.qza
03-repseqs/dada2-pe/aligned_dna_sequences.fasta
03-repseqs/dada2-pe/aligned_dna_sequences_gaps.txt
03-repseqs/dada2-pe/aligned_dna_sequences_gaps_describe.tsv
04-diversity/dada2-pe/taxa_barplot.qzv
04-diversity/dada2-pe/unweighted_unifrac_pcoa_results.qza
04-diversity/dada2-pe/observed_otus_vector.qza
04-diversity/dada2-pe/unweighted_unifrac_emperor.qzv
04-diversity/dada2-pe/weighted_unifrac_pcoa_results.qza
04-diversity/dada2-pe/weighted_unifrac_distance_matrix.qza
04-diversity/dada2-pe/unweighted_unifrac_distance_matrix.qza
04-diversity/dada2-pe/bray_curtis_pcoa_results.qza
04-diversity/dada2-pe/bray_curtis_distance_matrix.qza
04-diversity/dada2-pe/jaccard_pcoa_results.qza
04-diversity/dada2-pe/shannon_vector.qza
04-diversity/dada2-pe/jaccard_distance_matrix.qza
04-diversity/dada2-pe/rarefied_table.qza
04-diversity/dada2-pe/bray_curtis_emperor.qzv
04-diversity/dada2-pe/faith_pd_vector.qza
04-diversity/dada2-pe/jaccard_emperor.qzv
04-diversity/dada2-pe/evenness_vector.qza
04-diversity/dada2-pe/weighted_unifrac_emperor.qzv
04-diversity/dada2-pe/alpha_rarefaction.qzv
```

##### dada2_pe_stats (step 4.1)

```
04-diversity/dada2-pe/unweighted_unifrac_group_significance.qzv
```

##### dada2_pe_report (step 5)

```
01-imported/fastq_pe_count.csv
01-imported/fastq_pe_count_describe.tsv
05-reports/processing.txt
05-reports/report_dada2-pe.txt
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
* Were all my samples sequenced in the same sequencing run? (Rule check_illumina_run will check for this.)
* Do I have long enough sequening to do paired-end analysis, or do I have to do single-end analysis only?
* What sequence pre-processing has been done already: Demultiplexing? Quality filtering and FastQC? Primer removal? Merging of paired reads?

#### Assess sample set and metadata

* Is my metadata file complete? Are the relevant parameters of my dataset present as numeric or categorical variables?
* Do I have enough samples in each group of key metadata categories to determine an effect?

#### Format metadata and sequence data

* Is my metadata file properly formatted? See [Metadata in QIIME 2](https://docs.qiime2.org/2019.7/tutorials/metadata/), the [EMP Metadata Guide](http://www.earthmicrobiome.org/protocols-and-standards/metadata-guide/), and [QIIMP](https://qiita.ucsd.edu/iframe/?iframe=qiimp) for help formatting your metadata.
* Is my sequence data demultiplexed, in `.fastq.gz` format, and described in a QIIME 2 fastq manifest file? See [Fastq Manifest Formats](https://docs.qiime2.org/2019.7/tutorials/importing/#fastq-manifest-formats) from QIIME 2 for instructions for creating this file.
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

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
