<img src="https://upload.wikimedia.org/wikipedia/commons/0/00/Tourmaline-121240.jpg" height=200> <img src="http://melissabessmonroe.com/wp-content/uploads/2014/03/20140303_TourmalineSurfPark128.jpg" height=200>

# tourmaline

Amplicon sequencing is a metagenetics method whereby a single DNA locus in a community of organisms is PCR-amplified and sequenced. Tourmaline is an amplicon sequence processing workflow for Illumina sequence data that uses [QIIME 2](https://qiime2.org) and the software packages it wraps. Tourmaline manages commands, inputs, and outputs using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system.

Two methods of amplicon sequence processing are supported, both of which generate ASVs (amplicon sequence variants, which approximate the "true" or "exact" sequences in a sample) rather than OTUs (operational taxonomic units, which blur sequencing errors and microdiversity through clustering):

* [Deblur](https://github.com/biocore/deblur) is a greedy deconvolution algorithm based on known Illumina read error profiles ([Amir et al., 2017](https://doi.org/10.1128/mSystems.00191-16)).
* [DADA2](https://github.com/benjjneb/dada2) implements a quality-aware model of Illumina amplicon errors to infer sample composition by dividing amplicon reads into partitions consistent with the error model ([Callahan et al., 2016](https://doi.org/10.1038/nmeth.3869)).

Tourmaline is an alternative amplicon pipeline to [Banzai](https://github.com/jimmyodonnell/banzai), which was developed for [MBON](https://github.com/marinebon/MBON) and uses [Swarm](https://github.com/torognes/swarm) for OTU picking.

## Step 1. Data assessment

Consider your amplicon locus. Ask:

* What is the locus being amplified, and what are the primer sequences?
* How much sequence variation is expected for this locus and dataset?
* Is the expected sequence variation enough to answer my question?
* What is the expected amplicon size for this locus and dataset?

Evaluate the sequence data to determine the best parameters for processing. Ask:

* What type and length of sequencing was used? (e.g., MiSeq 2x150bp)
* Do I have long enough sequening to do paired-end analysis, or do I have to do single-end analysis only?
* What sequence pre-processing has been done already: Demultiplexing? Quality filtering and FastQC? Merging of paired reads?

Evaluate the sample set and sample metadata. Ask:

* Is my metadata file properly formatted? See the [QIIME 2 documentation](https://docs.qiime2.org/2018.6/tutorials/metadata/) and the metadata tool [QIIMP](https://qiita.ucsd.edu/iframe/?iframe=qiimp).
* Is my metadata file complete? Are the relevant parameters of my dataset present as numeric or categorical variables?
* Do I have enough samples in each group of key metadata categories to determine an effect?

## Step 2. Data preparation

Preprocess and format sequence data and metadata for QIIME 2 processing.

* Amplicon sequence data: preprocess, quality filter, and import into QIIME 2 artifact.
* Reference sequence data: format and import into QIIME 2 artifact.
* Metadata: format and import into QIIME 2 artifact.

## Step 3. Denoising (ASV picking)

Run Deblur or DADA2 to generate ASV feature tables and representative sequences. This is akin to OTU picking.

* Denoise amplicon data using Deblur or DADA2 to generate ASV feature tables (BIOM). Use paired-end mode (DADA2 only) if sequence and amplicon lengths permit.
* Rarefy tables to even depth.
* Summarize tables and representative sequences (ASVs).

## Step 4. ASV curation: phylogeny and taxonomy

Generate a phylogenetic tree of ASV sequences, and identify the taxonomy (phylum, class, order, etc.) of each ASV.

* Build a phylogenetic tree of ASV sequences, or insert ASV sequences into an existing tree for your amplicon locus.
* Assign taxonomy to ASVs using a reference database for your amplicon locus.

## Step 5. Core analyses: alpha/beta diversity and taxonomic profiles

* Alpha diversity: diversity metrics (evenness, Shannon, Faith's PD, observed sequences), alpha rarefaction, alpha group significance.
* Beta diversity: distance matrices (un/weighted UniFrac, Bray-Curtis, Jaccard), principal coordinates, Emperor plots, beta group significance.
* Taxonomy barplots.

## Step 6. Quality control

After completing processing and core analyses, determine if the results make sense. Ask:

* How many sequences did I start with, and how many are left after denoising?
* Are the representative sequences of similar length or of very different lengths?
* Do the sequence alignment and tree look reasonable?
* Do samples cluster in an expected way in PCoA space?
* Do the taxonomic profiles match expected taxonomic compositions?
