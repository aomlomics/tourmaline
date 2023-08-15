## Clone the Tourmaline directory

*Download the Tourmaline directory and files.*

Every time you run Tourmaline on a new dataset—or with new parameters on the same dataset—you will clone (download) a fresh copy of the Tourmaline GitHub repository (directory and contents) to your computer. This allows Snakemake to run smoothly without having to account for previous inputs and outputs. If you don't want to save the output from a previous run, it is fine to delete the output and re-run the Tourmaline commands in the same directory.

Navigate to your project directory and clone the Tourmaline repository there. If using the Docker container, a good place to use as your project directory is `/data`:

```
cd /data # or your project directory
git clone https://github.com/aomlomics/tourmaline.git
```

Rename the `tourmaline` directory to something else before working in it. We recommend adding the name of the project (e.g., `tourmaline-test` or `tourmaline-myproject`) and/or the date (e.g., `tourmaline-test-YYYYMMDD`). For example:

```
mv tourmaline tourmaline-test
```

## Initialize from an existing Tourmaline directory

*Setup a new Tourmaline run with the parameters from a previous run (optional).*

If you want to re-run Tourmaline on an existing dataset with a new set of parameters (and keep the output from your old parameters), you'll want to clone a fresh copy of the Tourmaline GitHub repository but then initialize it with all the parameters, files, and links you set up before. 

##### Script to initialize Tourmaline directory

After you clone a fresh copy of the Tourmaline GitHub repository, go inside that newly downloaded Tourmaline directory and run the script ***initialize_dir_from_existing_tourmaline_dir.sh***. This will copy `config.yaml` from an existing Tourmaline directory, remove the test files, then copy the data files and symlinks from the existing Tourmaline directory.

From the new Tourmaline top-level directory, run this:

```
scripts/initialize_dir_from_existing_tourmaline_dir.sh /path/to/tourmaline-existing
```

You may get error messages if some files don't exist, but it should have copied the files that were there. The files that will be copied from the existing directory to the new directory are:

```
config.yaml
00-data/manifest_pe.csv
00-data/manifest_se.csv
00-data/metadata.tsv
00-data/repseqs_to_filter_dada2-pe.tsv
00-data/repseqs_to_filter_dada2-se.tsv
00-data/repseqs_to_filter_deblur-se.tsv
01-imported/refseqs.qza
01-imported/reftax.qza
01-imported/classifier.qza
```

Check all the files. If everything looks good, you can breeze through the rest of the setup. You only need to make any new changes you want to the parameters in `config.yaml`. In some cases you may want to delete files you want to be regenerated.  

If you manually copy over output files from a previous Tourmaline run that you do not want to be regenerated (eg, 02-output-{method}-unfiltered), you should use the `cp -p` flag to preserve timestamps.

```
cp -rp tourmaline-old/02-output-dada2-pe-unfiltered/ tourmaline-new/
```

## Config file 

*Edit the config file to contain the file paths and parameters you want to use.*

The configuration file—named `config.yaml` and residing in the top-level directory—contains all of the parameters that are editable for a Tourmaline run. Always make sure the parameters are appropriate for your dataset. (Note: If you are comfortable editing Python code, it is possible to change more parameters by editing `Snakefile`.)

#### Test data: review the config file

The configuration file `config.yaml` comes ready to run with the test data. Note: Because the sequencing depth of the test dataset is very low (1000 sequences per sample) and the sequences are long (2x300 bp; 269 and 268 bp after trimming primers) relative to the sequenced amplicon length (~316 bp), several of the parameters in the config file may not be appropriate for a normal run of Tourmaline.

Below are the full contents of `config.yaml` with the default settings:

```
# config.yaml - configuration file for the Tourmaline Snakemake workflow
# Compatable with qiime2-2023.2
# User MUST edit these parameters before running their own data.
# Detailed instructions: https://github.com/aomlomics/tourmaline/wiki.

# METADATA FILE

# Metadata file must be named as follows (or use symbolic link):
# - 00-data/metadata.tsv

# Standardization is recommended through use of MIxS/MIMARKS (https://gensc.org/mixs/) column headers, for example:
# - submitted_to_insdc: {boolean}
# - investigation_type: [eukaryote|bacteria_archaea|plasmid|virus|organelle|metagenome|metatranscriptome|mimarks-survey|mimarks-specimen|misag|mimag|miuvig]
# - project_name: {text}
# - experimental_factor (text or EFO and/or OBI): {termLabel} {[termID]}|{text}
# - lat_lon (decimal degrees): {float} {float}

# FASTQ DATA (choose one option)

# Option 1 - if you want to start from manifest files and import per-sample fastq sequences,
# Name your manifest file(s) as follows (or use symbolic links):
# - paired-end sequences: 00-data/manifest_pe.csv
# - single-end sequences: 00-data/manifest_se.csv

# Option 2 - if your fastq sequences are already archived in qza format,
# Name your qza file(s) as follows (or use symbolic links):
# - paired-end sequences: 01-imported/fastq_pe.qza
# - single-end sequences: 01-imported/fastq_se.qza

# REFERENCE DATABASE (choose one option)

# Option 1 - if your reference database is not yet imported, 
# Name your reference fna and tsv files as follows (or use symbolic links):
# - reference sequences: 00-data/refseqs.fna
# - reference taxonomy: 00-data/reftax.tsv

# Option 2 - if your reference database is already archived in qza format,
# Name your your qza files as follows (or use symbolic links):
# - reference sequences: 01-imported/refseqs.qza
# - reference taxonomy: 01-imported/reftax.qza

# Provide a descriptive name for the reference database used, no spaces (e.g., original file name, source and version)
database_name: silva-138-99-515-806_q2-2021.2

# DENOISE

# DADA2 PAIRED-END
# For more info run: qiime dada2 denoise-paired --help

dada2pe_trunc_len_f: 240
dada2pe_trunc_len_r: 190
dada2pe_trim_left_f: 0
dada2pe_trim_left_r: 0
dada2pe_max_ee_f: 2
dada2pe_max_ee_r: 2
dada2pe_trunc_q: 2
dada2pe_pooling_method: independent
dada2pe_chimera_method: consensus
dada2pe_min_fold_parent_over_abundance: 1
dada2pe_n_reads_learn: 1000000
dada2pe_hashed_feature_ids: --p-hashed-feature-ids

# DADA2 SINGLE-END
# For more info run: qiime dada2 denoise-single --help

dada2se_trunc_len: 240
dada2se_trim_left: 0
dada2se_max_ee: 2
dada2se_trunc_q: 2
dada2se_pooling_method: independent
dada2se_chimera_method: consensus
dada2se_min_fold_parent_over_abundance: 1
dada2se_n_reads_learn: 1000000
dada2se_hashed_feature_ids: --p-hashed-feature-ids

# DEBLUR SINGLE-END
# For more info run: qiime deblur denoise-other --help

deblur_trim_length: 240
deblur_sample_stats: --p-sample-stats
deblur_mean_error: 0.005
deblur_indel_prob: 0.01
deblur_indel_max: 3
deblur_min_reads: 10
deblur_min_size: 2
deblur_hashed_feature_ids: --p-hashed-feature-ids

# TAXONOMIC CLASSIFICATION

# Taxonomic classification of the representative sequences.
# - method: choose from: naive-bayes, consensus-blast, consensus-vsearch
# - parameters: add arbitrary parameters or at least one parameter; see "qiime feature-classifier COMMAND --help"
# For more info run: qiime feature-classifier [fit-classifier-naive-bayes|classify-sklearn|classify-consensus-blast|classify-consensus-vsearch] --help
classify_method: consensus-vsearch
classify_parameters: --verbose

# Table of feature counts per sample collapsed by chosen taxonomic level.
# taxa_level - taxonomic level at which the features should be collapsed. Default 7
# For more info run: qiime taxa collapse --help
classify_taxalevel: 7

# MULTIPLE SEQUENCE ALIGNMENT

# Multiple sequence alignment of the representative sequences.
# - method: choose from: muscle v5, clustalo, mafft
# - muscle_iters: number of refinement iterations (integer, default 100)
# For more info run: muscle, clustalo --help, qiime alignment [mafft|mask] --help

alignment_method: muscle
alignment_muscle_iters: 50

# OUTLIER DETECTION

# Representative sequence outlier detection using odseq.
# - distance_metric: choose from: linear, affine
# - bootstrap_replicates: number of bootstrap replicates (integer)
# - threshold: probability to be at right of the bootstrap scores distribution (float)
# For more info see: https://www.bioconductor.org/packages/release/bioc/html/odseq.html

odseq_distance_metric: linear
odseq_bootstrap_replicates: 100
odseq_threshold: 0.025

# SUBSAMPLING (RAREFACTION)

# Subsample (rarefy) data to have an even number of observations per sample.
# - core_sampling_depth: subsampling depth for core diversity analyses
# - alpha_max_depth: maximum subsampling depth for alpha diversity rarefaction analysis

core_sampling_depth: 500
alpha_max_depth: 500

# BETA GROUP SIGNIFICANCE

# Statistical test for difference between samples grouped by a metadata variable (column).
# - column: metadata column to test beta group significance; this column must appear in the metadata file 00-data/metadata.tsv
# - method: choose from: permanova, anosim, permdisp
# - pairwise: choose from: --p-no-pairwise, --p-pairwise (can be very slow)

beta_group_column: region
beta_group_method: permanova
beta_group_pairwise: --p-pairwise

# DEICODE BETA DIVERSITY

# Robust Aitchison PCA (and biplot ordination) with automatically estimated underlying low-rank structure.
# Parameters are set to recommended defaults:
# - min_sample_count: minimum sum cutoff of samples across all features
# - min_feature_count: minimum sum cufoff of features across all samples
# - min_feature_frequency: minimum percentage of samples a feature must appear with value greater than zero
# - max_iterations: number of iterations to optimize the solution
# - num_features: number of most important features (arrows) to display in the biplot ordination
# For more info run: qiime deicode auto-rpca --help

deicode_min_sample_count: 500
deicode_min_feature_count: 10
deicode_min_feature_frequency: 0
deicode_max_iterations: 5
deicode_num_features: 5

# REPORT THEME

# Report theme for html version.
# Choose from: github, gothic, newsprint, night, pixyll, whitey.

report_theme: github

# FILTERING

# Filtering is implemented by filtering commands only, eg "snakemake dada2_pe_report_filtered", and is applied to each of
# representative sequences, taxonomy, and feature table.

# FILTER SAMPLES BY ID
# The file of sample IDs to be removed must have the sample IDs in the 1st column with a header line. Any number of additional columns can be added.
# The file must be named 00-data/samples_to_filter_{method}.tsv, where method is dada2-pe, dada2-se, or deblur-se.
# For more info run: qiime feature-table filter-samples --help
# TO SKIP FILTERING BY SAMPLE ID: provide a file with only headers (default) or don't run filtering commands.

# FILTER SAMPLES BY METADATA
# Samples can be removed or retained based on their metadata, using SQLite WHERE-clause syntax inside double quotes.
# For more info run: qiime feature-table filter-samples --help
# TO SKIP FILTERING BY SAMPLE METADATA: provide "none" (default) or don't run filtering commands.
# EXAMPLE: "[region]='Open Water'"
metadata_filter: none

# FILTER BY FEATURE ID
# The file of feature IDs to be filtered must have a header line with two columns: 1. "featureid", 2. anything.
# The file must be named 00-data/repseqs_to_filter_{method}.tsv, where method is dada2-pe, dada2-se, or deblur-se.
# If filtering outliers, just copy 02-output-{method}-{filter}/02-alignment-tree/repseqs_to_filter_outliers.tsv to the above filename.
# Add any additional feature IDs to be filtered (no duplicates allowed).
# For more info run: qiime feature-table [filter-seqs|filter-features] --help
# TO SKIP FILTERING BY FEATURE ID: provide only nonsense feature IDs in the above files (default) or don't run filtering commands.

# FILTER BY TAXONOMY
# Features with taxonomy containing these terms will be filtered.
# Separate terms with commas (e.g., mitochondria,chloroplast,eukaryota,unassigned).
# Terms are not case-sensitive.
# For more info run: qiime taxa filter-seqs --help
# TO SKIP FILTERING BY TAXONOMY: provide a nonsense term or don't run filtering commands.

exclude_terms: eukaryota,archaea,mitochondria,chloroplast,unassigned

# FILTER BY LENGTH
# Set minimum and maximum sequence lengths to filter representative sequences by.
# Limits are inclusive, ie, greater than or equal to minimum, less than or equal to maximum.
# For more info run: qiime feature-table filter-seqs --help
# TO SKIP FILTERING BY LENGTH: set values to extreme values, eg (0, 10000) or don't run filtering commands.

repseq_min_length: 0
repseq_max_length: 260

# FILTER BY ABUNDANCE & PREVALENCE
# Set minimum abundance/prevalence limits for filtering.
# Values are floats range(0,1)
# Limit is inclusive, ie, greater than or equal to minimum. Samples with frequency of 0 after filtering will also be removed.
# For more info run: qiime feature-table filter-features-conditionally --help
# TO SKIP FILTERING BY ABUNDANCE/PREVALENCE: set values to 0 or don't run filtering commands.

repseq_min_abundance: 0.01
repseq_min_prevalence: 0.1

# THREADS

# Max number of threads for individual rules.
# Threads used will be the lower of this and snakemake parameter --cores.
# Parameter other_threads is used for all other rules, regardless if they can use multiple threads, 
# because it prevents multiple rules from running simultaneously with --cores >1.

dada2pe_threads: 8
dada2se_threads: 8
deblur_threads: 8
alignment_threads: 8
feature_classifier_threads: 8
phylogeny_fasttree_threads: 8
diversity_core_metrics_phylogenetic_threads: 8
other_threads: 8
```

#### Your data: edit parameters customized to your data and analysis

Edit the configuration file `config.yaml` to change DADA2/Deblur parameters, subsampling (rarefaction) depth, multiple sequence alignment method, taxonomic classification method, and others parameters. You cannot change the input filenames in the sections "METADATA FILE", "FASTQ DATA" and "REFERENCE DATABASE"; you can either rename your files to match the standard filenames or use symbolic links to point your files to those filenames. You should consider changing the following parameters:

| Parameter                                                    | Description                                                  | Recommendation                                               | Help                                                         |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| `dada2pe_trunc_len_f` `dada2pe_trunc_len_r` `dada2se_trunc_len` | Truncate bases (integer) from the 3' (right) ends of reads in DADA2. | Choose values that maximize length but remove low-quality ends. Note that DADA2 paired-end mode requires a minimum overlap of 12 bp to merge Read 1 and Read 2. See the section below "Sequence quality control and choice of truncation length" for instructions on using the included script *fastqc_per_base_sequence_quality_dropoff.py*. | [dada2 denoise-paired](https://docs.qiime2.org/2023.2/plugins/available/dada2/denoise-paired/); [dada2 denoise-single](https://docs.qiime2.org/2023.2/plugins/available/dada2/denoise-single/) |
| `dada2pe_trim_left_f` `dada2pe_trim_left_f` `dada2se_trim_left` | Trim bases (integer) from the 5' (left) ends of reads in DADA2. | Depending on your amplicon sequencing method, and if trimming was not done prior to running Tourmaline, you may have primer sequences, indexes, and/or adapters on the 5' ends of your reads. If so, set this parameter to remove those bases. If not, set this parameter to zero. Note that 5' trimming (this parameter) is done *after* 3' truncation (above parameter). | [dada2 denoise-paired](https://docs.qiime2.org/2023.2/plugins/available/dada2/denoise-paired/); [dada2 denoise-single](https://docs.qiime2.org/2023.2/plugins/available/dada2/denoise-single/) |
| `deblur_trim_length`                                         | Truncate bases (integer) from the 3' (right) ends of reads in Deblur. | Choose values that maximize length but remove low-quality ends. See the section below "Sequence quality control and choice of truncation length" for instructions on using the included script *fastqc_per_base_sequence_quality_dropoff.py*. | [deblur denoise-other](https://docs.qiime2.org/2023.2/plugins/available/deblur/denoise-other/) |
| `dada2pe_pooling_method` `dada2se_pooling_method`            | DADA2 pooling method.                                        | Choose `pseudo` for pseudo-pooling or `independent` for no pooling. | [dada2 denoise-paired](https://docs.qiime2.org/2023.2/plugins/available/dada2/denoise-paired/); [dada2 denoise-single](https://docs.qiime2.org/2023.2/plugins/available/dada2/denoise-single/) |
| `dada2pe_chimera_method` `dada2se_chimera_method`            | DADA2 chimera method.                                        | Choose `pooled` if pseudo-pooling otherwise `consensus` or `none`. | [dada2 denoise-paired](https://docs.qiime2.org/2023.2/plugins/available/dada2/denoise-paired/); [dada2 denoise-single](https://docs.qiime2.org/2023.2/plugins/available/dada2/denoise-single/) |
| `alignment_method`                                           | Multiple sequence alignment method.                          | Choose `muscle` or `clustalo` for best accuracy or `mafft` for faster results. | [muscle](https://www.drive5.com/muscle/); [clustalo](http://www.clustal.org/omega/); [mafft](https://docs.qiime2.org/2023.2/plugins/available/alignment/) |
| `classify_method`                                            | Taxonomic classification method.                             | Choose `naive-bayes` for best accuracy or `consensus-blast` for faster results. | [feature-classifier](https://docs.qiime2.org/2023.2/plugins/available/feature-classifier/) |
| `exclude_terms`                                              | Filter terms (taxa) from taxonomy.                           | Specify terms (comma-separated, no spaces) to find in taxonomy and filter out (case-insensitive), or provide a nonsense term to skip this step when filtering. | [taxa filter-seqs](https://docs.qiime2.org/2023.2/plugins/available/taxa/filter-seqs/) |
| `repseq_min_length` `repseq_max_length`                      | Set minimum and maximum sequence lengths to filter representative sequences by. | Limits are inclusive, i.e., sequences will be retained if greater than or equal to minimum, less than or equal to maximum. Leave defaults (0, 10000) to do no filtering. | link                                                         |
| `repseq_min_abundance` `repseq_min_prevalence`               | set minimum abundance and prevalence limits to filter representative sequences by. | Limit is inclusive, i.e., sequences will be retained if greater than or equal to minimum. Leave default (0) to do no filtering. |                                                              |
| `odseq_distance_metric`                                      | Distance metric for odseq.                                   | Choose metric from: `linear`, `affine`.                      | [odseq](https://www.bioconductor.org/packages/release/bioc/html/odseq.html) |
| `odseq_bootstrap_replicates`                                 | Number (integer) of bootstrap replicates for odseq.          | Choose more replicates for more robust detection of outliers, fewer replicates for faster processing. | [odseq](https://www.bioconductor.org/packages/release/bioc/html/odseq.html) |
| `odseq_threshold`                                            | Threshold (float) for bootstrap probability distribution for odseq. | Probability to be at the right of the bootstrap scores distribution when computing outliers. Tune this parameter depending on the diversity and occurrence of outliers in the MSA. | [odseq](https://www.bioconductor.org/packages/release/bioc/html/odseq.html) |
| `core_sampling_depth`                                        | Rarefaction depth (integer) for core diversity metrics.      | Choose a value that balances sequencing depth (more is better) with number of samples retained (more is better). | [diversity core-metrics-phylogenetic](https://docs.qiime2.org/2023.2/plugins/available/diversity/core-metrics-phylogenetic/) |
| `alpha_max_depth`                                            | Rarefaction depth (integer) for alpha rarefaction.           | Choose a value that balances sequencing depth (more is better) with number of samples retained (more is better). | [diversity alpha-rarefaction](https://docs.qiime2.org/2023.2/plugins/available/diversity/alpha-rarefaction/) |
| `beta_group_column`                                          | Column (text) in your metadata to test beta-diversity group significance. | Choose a category that may differentiate your samples. This analysis can be rerun with different columns by renaming the output file and changing the value in *config.yaml* before running again. | [diversity beta-group-significance](https://docs.qiime2.org/2021.2/plugins/available/diversity/beta-group-significance/) |
| `report_theme`                                               | HTML report theme.                                           | Choose from: `github`, `gothic`, `newsprint`, `night`, `pixyll`, `whitey`. | [Typora theme gallery](https://theme.typora.io/)             |

## Reference database

*Format reference sequences and taxonomy for QIIME 2, if not done already.*

Reference sequences and taxonomy are used by QIIME 2 to assign taxonomy to the representative sequences you get from DADA2 and Deblur. Any marker gene (amplicon locus) can be used with Tourmaline. The only requirement is that the reference sequence (.fasta) and taxonomy (.tsv) files are formatted for QIIME 2: the sequence IDs in the two files must match, and the taxonomy must have two tab-delimited columns where the second column is taxonomic levels delimited by semicolons. Here is what they should look like (shortened for simplicity):

`refseqs.fna`

```
>AB302407.1.2962
TCCGGTTGATCCTGCCGGACCCGACCGCTATCGGGGTAGGGCTAAGCCATGCGA
>KU725476.45629.48552
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCATGCTTAACACATGCAAG
```

`reftax.tsv`

```
AF352532.1.1501	D_0__Bacteria;D_1__Patescibacteria;D_2__Parcubacteria;D_3__GW2011-GWA2-46-7
KM650414.1.1439	D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia;D_3__Clostridiales
```

The reference sequence (.fasta) and taxonomy (.tsv) files, and their derived QIIME 2 archives (.qza), can be used for multiple studies and need only exist in one place on your computer. Therefore, it is best to store these files in a special directory for databases—for example, in a directory named `~/databases/qiime2`—then create symbolic links to them in the respective project subdirectories (you'll see how to do this below).

When you run Snakemake, if you get an error like, "Missing input files for rule import_ref_seqs: 00-data/refseqs.fna", you can fix this either by placing the reference sequence (.fasta) and taxonomy (.tsv) files in `00-data` or by placing their derived QIIME 2 archives (.qza) in `01-imported`. You will do the latter below.

#### Test data: download reference database and create symbolic links

The small test dataset included with Tourmaline is a 16S rRNA amplicon dataset sequenced using the 515F-806R primers, so we will download a 16S+18S reference database trimmed to this region (Silva 138 SSURef NR99 515F/806R region sequences and taxonomy). Please note the citations and license for this dataset [here](https://docs.qiime2.org/2023.2/data-resources/#silva-16s-18s-rrna).

```
cd ~/databases/qiime2
wget https://data.qiime2.org/2021.2/common/silva-138-99-seqs-515-806.qza
wget https://data.qiime2.org/2021.2/common/silva-138-99-tax-515-806.qza
```

Now navigate to the directory `01-imported` inside your new Tourmaline directory:

```
cd /path/to/tourmaline-test/01-imported
```

Create symbolic links to the database files:

```
ln -s ~/databases/qiime2/silva-138-99-seqs-515-806.qza refseqs.qza
ln -s ~/databases/qiime2/silva-138-99-tax-515-806.qza reftax.qza
```

#### Your data: format and import your reference database

There are three options for setting up your QIIME 2 reference database files. In each case, the underlying reference sequence (.fna) and taxonomy (.tsv) text files have to follow the QIIME 2 format.

A QIIME 2 reference sequence (.fna) file looks like this:

```
>AB302407
TCCGGTTGATCCTGCCGGACCCGACCGCTATCGGGGTAGGGCTAAGCCATGCGAGTCGCGCGCCCGGGGGCGCCCGGGAG
>KU725476
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCATGCTTAACACATGCAAGTCGGACGTTGTCTTCAAATTAGAATA
>KP109803
AGAGTTTGATCCTGGCTCAGAATGAACGCTAGCGATATGCTTAACACATGCAAGTCGAACGTTGTCTTCAAATTAGAATA
```

A QIIME 2 reference taxonomy (.tsv) file looks like this (sequence IDs must match those in the sequence file, not necessarily in the same order):

```
AB302407<TAB>D_0__Archaea;D_1__Crenarchaeota;D_2__Thermoprotei;D_3__Thermoproteales;D_4__Thermoproteaceae;D_5__Pyrobaculum;D_6__Pyrobaculum sp. M0H
KU725476<TAB>D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rickettsiales;D_4__Mitochondria;D_5__Sphagnum girgensohnii;D_6__Sphagnum girgensohnii
KP109803<TAB>D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rickettsiales;D_4__Mitochondria;D_5__Sphagnum palustre;D_6__Sphagnum palustre
```

##### Option 1 (import reference database using Tourmaline)

Start with the reference sequence and taxonomy files and let Tourmaline create the QIIME 2 archives (.qza). Place the text files (or symbolic links to them) in `00-data` and Tourmaline will create the latter for you in `01-imported`. For example, if your Tourmaline directory is `/path/to/tourmaline-myproject` and your (QIIME 2-ready) reference sequences and taxonomy are in `~/databases/qiime2/16s`, you could create symbolic links them using the commands below (you can change the filenames `refseqs.fna` and `reftax.tsv` in `config.yaml`, but we don't recommend it):

```
cd /path/to/tourmaline-myproject/00-data
ln -s ~/databases/qiime2/16s/16s_refseqs.fna refseqs.fna
ln -s ~/databases/qiime2/16s/16s_reftax.tsv reftax.tsv
```

##### Option 2 (import reference database on your own or use pre-imported database)

Import your reference database files to .qza files yourself before running Tourmaline, or use the .qza files from a previous run or the internet (e.g., the Silva .qza files mentioned above).

To import your own reference database files, use these commands:

```
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path 00-data/refseqs.fna \
--output-path 01-imported/refseqs.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path 00-data/reftax.tsv \
--output-path 01-imported/reftax.qza
```

If you already have the .qza files, you can put them in your databases directory (somewhere on your computer) and create symbolic links to them in `01-imported`. Tourmaline will see these files and skip the step that creates them (you can change the filenames `refseqs.qza` and `reftax.qza` in `config.yaml`, but we don't recommend it):

```
cd /path/to/tourmaline-myproject/01-imported
ln -s ~/databases/qiime2/16s/16s_refseqs.qza refseqs.qza
ln -s ~/databases/qiime2/16s/16s_reftax.qza reftax.qza
```

**Advanced:** If you have extracted a specific region or created a naive Bayes classifier in a previous run, you can place these in `01-imported` also:

```
ln -s ~/databases/qiime2/16s/16s_refseqs_515f_806r.qza refseqs_extracted.qza
ln -s ~/databases/qiime2/16s/16s_classifier.qza classifier.qza
```

For Snakemake to work with these symbolic links, you may have to run `snakemake --cleanup-metadata <filenames>` on them first.

## Sequence files

*List fastq.gz files in a QIIME 2 FASTQ manifest tsv file(s), or use your pre-imported FASTQ QIIME 2 archive file(s).* 

Tourmaline requires that your amplicon sequence data is already demultiplexed. We strongly recommend checking if your sequences have primers prior to running Tourmaline, and either [trimming them yourself](https://github.com/aomlomics/tourmaline/wiki/3-Setup#remove-primers-if-necessary) or choosing the appropriate truncation and trimming parameters in the config file.  

There are two options to provide your sequence data:

Option 1 is to keep your data as per-sample fastq.gz files with a FASTQ manifest file and let Tourmaline create the FASTQ archive (.qza) file(s). The FASTQ manifest files—by default named `manifest_se.csv` and `manifest_pe.csv` and found in `00-data`—tell QIIME 2 where to find your FASTQ files and what samples they correspond to. See [QIIME 2 Docs](https://docs.qiime2.org/2023.2/tutorials/importing/#fastq-manifest-formats) for more information on FASTQ manifest files. Hint: Gzipped files like fastq.gz files can be viewed using `zcat` (or `gzcat` on Mac) and `zless`.

Option 2 is to use a pre-imported .qza file (QIIME 2 archive), which already contains your demultiplexed sequences and a manifest file. 

#### Test data: edit the manifest file to match your filepath

We'll use Option 1 (FASTQ manifest and individual fastq.gz files) for the test data. The test sequence data were subsampled to 1000 sequences per sample and are provided in `00-data`. They are already demultiplexed and in ".fastq.gz" format, and the file names match those in the FASTQ manifest files and metadata file. To run the test data, unless you are using the Docker container, you must edit the manifest files to point to the absolute filepaths of the sequences in your local copy of `tourmaline` (which you renamed to `tourmaline-test`). For example, if the filepath of your project is `/Users/me/myproject`, these commands will set up the manifest files appropriately:

```
cd /Users/me/myproject/tourmaline-test/00-data
cat manifest_pe.csv | sed 's|/data/tourmaline|/Users/me/myproject/tourmaline-test|' > temp
mv temp manifest_pe.csv 
cat manifest_pe.csv | grep -v "reverse" > manifest_se.csv
```

If you are using the Docker container and you cloned `tourmaline` into `/data`, you don't need to modify the manifest files!

#### Your data: choose manifest+sequences or pre-imported sequence artifact

First, delete the test fastq.gz files that came in directory `00-data`.

##### Demultiplex (if necessary)

If your sequences are not already demultiplexed, demultiplex them, then gzip the resulting per-sample FASTQ files. See the [QIIME 2 Docs](https://docs.qiime2.org/2023.2/tutorials/importing/) for more information on importing and demultiplexing FASTQ files. If your sequences are in multiplexed EMP format, you can use the commands `qiime tools import` to import them as type `EMPPairedEndSequences`, then `qiime demux emp-paired` to demuliplex them.

##### Remove primers (if necessary)

If primers sequences are present in your sequences, you can use QIIME 2's Cutadapt command to trim the primers. The commands below show how to trim **16S V4-V5 primers** from paired-end data and then format the output for Tourmaline processing:

1. Load paired-end 16S FASTQ files into QIIME 2. FASTQ files are stored in a subfolder (paired_reads). This will result in a .qza file with information on 16S sequences.

    ```
    qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path paired_reads --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end-16S.qza
    ```

2. Run Cutadapt on the demultiplexed 16S reads to remove V4-V5 primers (change to your specific primers). 

    ```
    qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end-16S.qza --p-front-f GTGYCAGCMGCCGCGGTAA --p-adapter-f AAACTYAAAKRAATTGRCGG --p-front-r CCGYCAATTYMTTTRAGTTT --p-adapter-r TTACCGCGGCKGCTGRCAC --verbose --p-error-rate 0.2 --o-trimmed-sequences demux-trimmed-16S.qza
    ```

3. Export trimmed reads back to FASTQ file format in a new folder (16S-trimmed-reads). These reads can now be loaded into Tourmaline. And off you go with the protocol for the trimmed reads. 

    ```
    qiime tools export --input-path demux-trimmed-16S.qza --output-path 16S-trimmed-reads
    ```

If your amplicon length is shorter than your sequencing read length, you may need to do cutadapt in two steps, for example this strategy for 18S V9 primers with 250 bp paired-end sequencing:  

1. Load paired-end 18S FASTQ files into QIIME 2. FASTQ files are stored in a subfolder (paired_reads). This will result in a .qza file with information on 16S sequences.

    ```
    qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path paired_reads --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end-18S.qza
    ```

2. Run Cutadapt on the demultiplexed 18S reads to remove V9 primers from the ends of reads (change to your specific primers). 

    ```
    qiime cutadapt trim-paired \
    --i-demultiplexed-sequences demux-paired-end-18S.qza \
    --p-adapter-f AGTAGGTGAACCTGCAGAAGGATC #reverse comp of reverse primer\
    --p-adapter-r GACGGGCGGTGTGTAC #reverse comp of fwd primer\
    --p-match-read-wildcards \
    --p-match-adapter-wildcards \
    --verbose \
    --o-trimmed-sequences trimmed_remove_primers_wild.qza
    ```

3. Run cutadapt again on your trimmed reads, this time trimming primers off the front and only keeping trimmed sequences containing primers.

    ```
    qiime cutadapt trim-paired \
    --i-demultiplexed-sequences trimmed_remove_primers_wild.qza \
    --p-front-f GTACACACCGCCCGTC #forward primer\
    --p-front-r TGATCCTTCTGCAGGTTCACCTAC #reverse primer\
    --p-match-read-wildcards \
    --p-match-adapter-wildcards \
    --p-discard-untrimmed \
    --verbose \
    --o-trimmed-sequences trimmed_remove_primers_wild_2.qza
    ```

4. Export trimmed reads back to FASTQ file format in a new folder (18S-trimmed-reads). These reads can now be loaded into Tourmaline. And off you go with the protocol for the trimmed reads. 

    ```
    qiime tools export --input-path trimmed_remove_primers_wild_2.qza --output-path 18S-trimmed-reads
    ```

##### Option 1 (FASTQ manifest and individual fastq.gz files)

Find the location of your sequence data on your computer or external drive. It is recommended that you leave the sequences where they are and just tell Tourmaline where to find them using the FASTQ manifest file(s) in the next step.

Using the sample names in your metadata file and the absolute filepaths of the forward and reverse demultiplexed gzipped sequence files for each sample, create your FASTQ manifest file. Ensure the sample IDs match those in your metadata file. The included script `create_manifest_from_metadata.py` (see below) can help with this.

Tip: Create and finalize your paired-end manifest (*manifest_pe.csv*) first, then run this command to generate the single-end manifest (*manifest_se.csv*):

```
cat manifest_pe.csv | grep -v "reverse" > manifest_se.csv
```

Two scripts are provided to help you create your FASTQ manifest files:

***create_manifest_from_fastq_directory.py*** – Create a FASTQ manifest file from a directory of FASTQ files. This script makes the following assumptions: the first characters after the sample names are "\_S[0-9]{1,3}""; only Lane 1 data is present ("\_L001"); R1 and R2 files are both present; and only one file (with suffix "\_001") is present for each of R1 and R2.

```
scripts/create_manifest_from_fastq_directory.py /path/to/fastq \
00-data/manifest_pe.csv 00-data/manifest_se.csv
```

***match_manifest_to_metadata.py*** – Given a metadata file and a FASTQ manifest file, generate two new manifest files (paired-end and single-end) corresponding to the samples in the metadata file. This can be useful if you demultiplexed from "Earth Microbiome Project (EMP) protocol" format for paired-end reads and unzipped the output to get the manifest file.

```
scripts/match_manifest_to_metadata.py 00-data/metadata.tsv /path/to/demux/data/MANIFEST \
00-data/manifest_pe.csv 00-data/manifest_se.csv
```

##### Option 2 (pre-imported sequences as .qza file)

Put the .qza files in `01-imported`. It is recommended that you name them (or use symbolic links to name them) `fastq_pe.qza` (paired-end data) and `fastq_se.qza` (single-end data), but you can change these filenames in `config.yaml`.

If your .qza file is in the correct format, the output of the command `qiime tools peek FILE` should be:

Paired-end data:

```
Type:        SampleData[PairedEndSequencesWithQuality]
Data format: SingleLanePerSamplePairedEndFastqDirFmt
```

Single-end data:

```
Type:        SampleData[SequencesWithQuality]
Data format: SingleLanePerSampleSingleEndFastqDirFmt
```

## Sequence QC and truncation length

After your sequences are demultiplexed to per-sample .fastq.gz files, you can (optionally) run some diagnostics on them:

* FastQC to get error profiles and lots of other information about your sequences.
* MultiQC to collate and summarize the FastQC results.
* Helper script *fastqc_per_base_sequence_quality_dropoff.py* to help choose appropriate truncation (DADA2) and trim (Deblur) lengths for your configuration file.

##### FastQC & MultiQC

Installation note: This section requires a separate Conda environment, `multiqc`, which you can create and activate with these commands:

```
conda create -n multiqc -c bioconda fastqc multiqc
conda activate multiqc
```

First, create directories for FastQC and MultiQC output for Read 1 and Read 2:

```
cd /path/to/fastq
mkdir fastqc-R1
mkdir fastqc-R2
```

Then run FastQC on the forward and the reverse reads (change the .fastq.gz filename wildcards as necessary):

```
fastqc *_R1_001.fastq.gz -o fastqc-R1
fastqc *_R2_001.fastq.gz -o fastqc-R2
```

Then run MultiQC on the output of each FastQC run:

```
cd fastqc-R1
multiqc --export .
cd ../fastqc-R2
multiqc --export .
```

##### Choose sequence truncation lengths

Finally, reactivate your `qiime2-2023.2` environment and run ***fastqc_per_base_sequence_quality_dropoff.py***, which will determine the position where median per-base sequence quality drops below some fraction (default: 0.90, here: 0.85) of its maximum value. The output of this script can be used as your parameter setting for DADA2 truncation or Deblur trim values (use output from Read 1 for `dada2pe_trunc_len_f`, `dada2se_trunc_len`, and `deblur_trim_length`; use output from Read 2 for `dada2pe_trunc_len_r`):

```
cd /path/to/tourmaline
scripts/fastqc_per_base_sequence_quality_dropoff.py \
--input 00-data/fastq/fastqc-R1/multiqc_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt \
--cutoff 0.85
```

Note that DADA2 paired-end mode requires a minimum overlap of 12 bp to merge Read 1 and Read 2.

## Metadata file

*Format metadata file for QIIME 2.*

The metadata file—by default named `metadata.tsv` and found in `00-data`—is a tab-delimited text file (e.g., exported from Excel) with samples as rows and metadata categories as columns. It is also known as the mapping file.

#### Test data: metadata is ready to go

The metadata file that comes with Tourmaline is ready to go with the test data. It is complete, curated, and its sample names match the sample names in the FASTQ manifest file.

#### Your data: standardize and curate your metadata with meaningful categories

It's important to have rich and complete sample metadata before you begin your analyses. Your sample metadata should include basic collection information like `collection_timestamp`, `latitude`, and `longitude`, plus categories that describe treatment groups and environmental metadata relevant to your samples.

Sample preparation information is also part of your metadata and includes information about how the samples were sequenced.

The following are some suggested columns to include in your metadata for each project:

* `project_name`
* `experiment_design_description`
* `target_gene`
* `target_subfragment`
* `pcr_primers`
* `pcr_primer_names`
* `platform`
* `instrument_model`
* `run_center`
* `run_date`

The above columns follow the standards set by [Qiita](https://qiita.ucsd.edu/static/doc/html/gettingstartedguide/index.html). For additional help formatting your metadata see [Metadata in QIIME 2](https://docs.qiime2.org/2023.2/tutorials/metadata/), the [EMP Metadata Guide](http://www.earthmicrobiome.org/protocols-and-standards/metadata-guide/), [QIIMP](https://qiita.ucsd.edu/iframe/?iframe=qiimp), and NMDC's [Introduction to Metadata and Ontologies](https://microbiomedata.org/introduction-to-metadata-and-ontologies/).

##### Symlinks to metadata and FASTQ manifest files

It's a good idea to keep current versions of your metadata file (and perhaps your FASTQ manifest files) in a dedicated metadata directory in your project directory. Then make a symbolic link to these files in the `00-data` directory.  Use the code below to remove the test files and create links to your own metadata and FASTQ manifest files (the files in `/path/to/metadata` can be named however you like):

```
cd /path/to/tourmaline/00-data
rm metadata.tsv manifest_se.csv manifest_pe.csv
ln -s /path/to/metadata/metadata.tsv metadata.tsv
ln -s /path/to/metadata/manifest_se.csv manifest_se.csv
ln -s /path/to/metadata/manifest_pe.csv manifest_pe.csv
```

## Appendix: Helper scripts

As we have seen already, Tourmaline comes with helper scripts in the `scripts` directory. 

#### Setup

These scripts might come in handy in setting up your workflow.

***initialize_dir_from_existing_tourmaline_dir.sh*** – From the main directory of a newly cloned tourmaline directory, run this script to copy the config.yaml and Snakefile from an existing tourmaline directory, remove the test files, then copy the data files and symlinks from the existing tourmaline directory. Very useful when re-running an analysis on the same dataset. Simply clone a new copy of Tourmaline, then run this script to copy everything from the old Tourmaline directory to the new one, then make your desired changes to the parameters.

```
scripts/initialize_dir_from_existing_tourmaline_dir.sh /path/to/existing/tourmaline
```

***create_manifest_from_fastq_directory.py*** – Create a FASTQ manifest file from a directory of FASTQ files. This script makes the following assumptions: the first characters after the sample names are "\_S[0-9]{1,3}"; only Lane 1 data is present ("\_L001"); R1 and R2 files are both present; and only one file (with suffix "\_001") is present for each of R1 and R2.

```
scripts/create_manifest_from_fastq_directory.py /path/to/fastq/dir \
00-data/manifest_pe.csv 00-data/manifest_se.csv
```

***match_manifest_to_metadata.py*** – Given a metadata file and a FASTQ manifest file, generate two new manifest files (paired-end and single-end) corresponding to the samples in the metadata file. This can be useful if you demultiplexed from "Earth Microbiome Project (EMP) protocol" format for paired-end reads and unzipped the output to get the manifest file.

```
scripts/match_manifest_to_metadata.py 00-data/metadata.tsv /path/to/demux/data/MANIFEST \
00-data/manifest_pe.csv 00-data/manifest_se.csv
```

***fastqc_per_base_sequence_quality_dropoff.py*** – Determine the position where median per base sequence quality drops below some fraction (default: 0.90) of its maximum value. This is useful for defining 3' truncation positions in DADA2 (trunc-len) and Deblur (trim-length). This script should be run separately for Read 1 and Read 2 fastqc/multiqc output.

```
scripts/fastqc_per_base_sequence_quality_dropoff.py \
--input /path/to/fastqc-R1/multiqc_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt \
--cutoff 0.85
```

#### Wrapped

These scripts are wrapped by Tourmaline (used by the workflow), but some of them can be run separately.

***detect_amplicon_locus.py*** – Try to guess the amplicon locus of a fasta file based on the first four nucleotides. Used by rule `repseq_detect_amplicon_locus`.

```
scripts/detect_amplicon_locus.py \
-i 02-output-dada2-pe/00-table-repseqs/repseqs.fasta \
> 02-output-dada2-pe/00-table-repseqs/repseqs_amplicon_type.txt
```

***fastaLengths.pl*** – Calculate the length of every sequence in a fasta file. Used by rule `repseq_length_distribution`.

```
scripts/fastaLengths.pl 02-output-dada2-pe/00-table-repseqs/repseqs.fasta \
> 02-output-dada2-pe/00-table-repseqs/repseqs_lengths.txt
```

***run_odseq.R*** – Determine which representative sequences are likely to be outliers using the R package *odseq*. Used by rule `alignment_detect_outliers`.

```
Rscript --vanilla scripts/run_odseq.R 02-output-dada2-pe/01-alignment-tree-taxonomy/aligned_dna_sequences.fasta affine 1000 0.025 02-output-dada2-pe/01-alignment-tree-taxonomy/aligned_dna_sequences_outliers.csv
```