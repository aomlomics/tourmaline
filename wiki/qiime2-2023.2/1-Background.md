## Methods

Tourmaline builds on the latest methods for analysis of microbial and eDNA amplicon sequence data. This section describes those methods and provides tutorials for some of them.

### Amplicon sequence variants

Amplicon sequencing (metabarcoding) is a method whereby a single DNA locus in a community of organisms is PCR-amplified and sequenced. Two methods of amplicon sequence processing are supported, both of which generate ASVs (amplicon sequence variants), which approximate the "true" or "exact" sequences in a sample, rather than OTUs (operational taxonomic units), which blur sequencing errors and microdiversity through clustering:

* [DADA2](https://github.com/benjjneb/dada2) implements a quality-aware model of Illumina amplicon errors to infer sample composition by dividing amplicon reads into partitions consistent with the error model ([Callahan et al., 2016](https://doi.org/10.1038/nmeth.3869)).
* [Deblur](https://github.com/biocore/deblur) is a greedy deconvolution algorithm based on known Illumina read error profiles ([Amir et al., 2017](https://doi.org/10.1128/mSystems.00191-16)).

### QIIME 2 

[QIIME 2](https://qiime2.org/) ([Bolyen et al., 2019](https://doi.org/10.1038/s41587-019-0209-9)) is one of the most popular amplicon sequence analysis software tools available. It supports both DADA2 and Deblur denoising algorithms as well as a variety of downstream diversity and statistical analyses and visualizations. [Click here for a tutorial on QIIME 2.](https://github.com/aomlomics/tutorials/tree/master/qiime2)

### Snakemake

[Snakemake](https://snakemake.readthedocs.io/en/stable/) is a workflow management software that allows for reproducible and scalable workflows for bioinformatics and other data analyses. It keeps track of input and output files, storing output files in a logical directory structure. It uses rules to define commands and only runs rules when they are required to produce the desired output. [Click here for a tutorial on Snakemake.](https://github.com/aomlomics/tutorials/tree/master/snakemake)

## Assessing your data

Before starting any bioinformatics workflow, it is important to assess your data and metadata to decide how they need to be formatted and inform your parameter choices. The questions below can help determine the best parameters for processing and to be able to evaluate the success of your completed workflow.

### Amplicon locus

* What is the locus being amplified, and what are the primer sequences?
* How much sequence variation is expected for this locus (and primer sites) and dataset?
* Is the expected sequence variation enough to answer my question?
* What is the expected amplicon size for this locus and dataset?

### Sequence data

* What type and length of sequencing was used? (e.g., MiSeq 2x150bp)
* Were all my samples sequenced in the same sequencing run? (Rule `check_illumina_run` will check for this.)
* Do I have long enough sequening to do paired-end analysis, or do I have to do single-end analysis only?
* Has the relevant sequence pre-processing been done already: Demultiplexing? Quality filtering and/or trimming? Primer removal? FastQC/MultiQC profiling before and/or after filtering/trimming? (Note on quality trimming: If you plan on using DADA2 for denoising, the developers recommend no quality trimming be done before running, as it confuses the error profiling algorithm of DADA2. Quality filtering including complete removal of erroneous sequences is still advisable.)

### Sample set and metadata

* Is my metadata file complete? Are the relevant parameters of my dataset present as numeric or categorical variables?
* Do I have enough samples in each group of key metadata categories to determine an effect?

## Workflow overview

The table below describes the basic steps of the Tourmaline workflow. Further instructions are provided in the sections [Setup](https://github.com/aomlomics/tourmaline/wiki/3-Setup) and [Run](https://github.com/aomlomics/tourmaline/wiki/4-Run).

In the file paths below, `{method}` is one of:

* `dada2-pe` (paired-end DADA2)
* `dada2-se` (single-end DADA2)
* `deblur-se` (single-end Deblur)

and `{filter}` is one of:

* `unfiltered` (representative sequences and feature table *are not* filtered by taxonomy or feature ID)
* `filtered` (representative sequences and feature table *are* filtered by taxonomy or feature ID)

| Step                                | Command                        | Output                                                       |
| ----------------------------------- | ------------------------------ | ------------------------------------------------------------ |
| Format input and configuration file | (ad hoc)                       | `config.yaml`, `Snakefile`, `00-data/metadata.tsv`, `00-data/manifest_se.tsv`, `00-data/manifest_pe.tsv`, `00-data/refseqs.fasta` or `01-imported/refseqs.qza`, `00-data/reftax.tsv` or `01-imported/reftax.qza` |
| Import data                         | `snakemake {method}_denoise`   | `01-imported/` (multiple files)                              |
| Denoising                           | `snakemake {method}_denoise`   | `02-output-{method}/00-table-repseqs/` (multiple files)      |
| Taxonomic assignment    | `snakemake {method}_taxonomy_{filter}` | `02-output-{method}/01-taxonomy` (multiple files) |
| Representative sequence curation    | `snakemake {method}_diversity_{filter}` | `02-output-{method}/02-alignment-tree` (multiple files) |
| Core diversity analyses             | `snakemake {method}_diversity_{filter}` | `02-output-{method}-{filter}/03-alpha-diversity/` `02-output-{method}-{filter}/04-beta-diversity/` (multiple files) |
| Report                              | `snakemake {method}_report_{filter}`    | `03-reports/report_{method}_{filter}.html`                            |

## Contact us 

* Questions? Join the conversation on [gitter](https://gitter.im/aomlomics/tourmaline).
* Have a feature request? Raise an issue on [GitHub](https://github.com/aomlomics/tourmaline/issues).

 