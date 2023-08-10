## Run Snakemake

Now we are ready to start running Tourmaline commands using the `snakemake` command. In the file paths below, `{method}` is one of: `dada2-pe` (paired-end DADA2), `dada2-se` (single-end DADA2), or `deblur-se` (single-end Deblur); and `{filter}` is one of `unfiltered` (no filtering of representative sequences and feature tables) or `filtered` (filtering of representative sequences and feature tables by taxonomy and/or feature ID).

Note that any of the commands below can be run with various options, including `--printshellcmds` to see the shell commands being executed and `--dryrun` to display which rules would be run but not execute them.

Tourmaline will automatically check if the required input files and parameters are provided properly before running, and will generate an error with feedback if they are not.

### Activate snakemake environment

All of the commands below will be run with your 'snakemake' environment activated. The snakemake commands will then activate the 'qiime2-2023.5' environment automatically as they are run.  

```
conda activate snakemake
```

#### Number of compute cores

You must provide the maximum number of threads you want with the `snakemake --cores` parameter. If the `config.yaml` file specifies more threads than the number specified by `--cores`, it will scale down to match the `--cores` value.

Example, running with 5 cores:

```
snakemake --use-conda dada2_pe_denoise --cores 5
```

Every snakemake command below must be run with the '--use-conda' option. 

### Unfiltered mode

The first command is `snakemake dada2_pe_denoise`, which imports the FASTQ files and the reference database (if not already present in directory `01-imported`), summarizes the FASTQ data, runs denoising using DADA2, and summarizes and visualizes the output feature table and representative sequences (for DADA2, the representative sequences are amplicon sequence variants or ASVs). 

```
snakemake --use-conda dada2_pe_denoise --cores 5
# or dada2_se_denoise or deblur_se_denoise
```

At this point, the user should determine if the DADA2 parameters need to be modified (and if so, delete the output files and rerun the *denoise* step) and when satisfied choose appropriate rarefaction depths for the parameters “alpha_max_depth” and “core_sampling_depth” in `config.yaml`. 

The second command is `snakemake dada2_pe_taxonomy_unfiltered`, which assigns taxonomy to the representative sequences using a Naive Bayes classifier or consensus BLAST method, visualizes the taxonomy, and generates an interactive taxa barplot. 

```
snakemake --use-conda dada2_pe_taxonomy_unfiltered --cores 5
# or dada2_se_taxonomy_unfiltered or deblur_se_taxonomy_unfiltered
```

The third command is `snakemake dada2_pe_diversity_unfiltered`, which will align and build a phylogenetic tree of the representative sequences, identify representative sequences that have unassigned taxonomy or are potential outliers, summarize and plot the representative sequence properties, perform alpha-rarefaction, and run alpha-diversity and beta-diversity analyses and group significance tests using a full suite of metrics. 

```
snakemake --use-conda dada2_pe_diversity_unfiltered --cores 5
# or dada2_se_diversity_unfiltered or deblur_se_diversity_unfiltered
```

The fourth and final command is `snakemake dada2_pe_report_unfiltered`, which will create a comprehensive HTML report of parameters, metadata, inputs, outputs, and visualizations in a single file. The report includes hyperlinks to QIIME 2 visualization files, which can be downloaded and drag-and-dropped into [view.qiime2.org](https://view.qiime2.org) for viewing. Filtering of representative sequences is provided by the *filtered* mode.

```
snakemake --use-conda dada2_pe_report_unfiltered --cores 5
# or dada2_se_report_unfiltered or deblur_se_report_unfiltered
```

### Filtered mode

After viewing the *unfiltered* results—the taxonomy summary and taxa barplot, the representative sequence summary plot and table, or the list of unassigned and potential outlier representative sequences—the user may wish to filter (remove) certain taxonomic groups or representative sequences. If so, the user should first check the following parameters and/or files:

* copy `2-output-dada2-pe-unfiltered/02-alignment-tree/repseqs_to_filter_outliers.tsv` to `00-data/repseqs_to_filter_dada2-pe.tsv` to filter outliers, or manually include feature IDs in `00-data/repseqs_to_filter_dada2-pe.tsv` to filter those feature IDs (change "dada2-pe" to "dada2-se" or "deblur-se" as appropriate);
*  `exclude_terms` in `config.yaml` – add taxa to exclude from representative sequences, if desired;
*  `repseq_min_length` and `repseq_max_length` in `config.yaml` – set minimum and/or maximum lengths for filtering representative sequences, if desired;
*  `repseq_min_abundance` and `repseq_min_prevalence` in `config.yaml` – set minimum abundance and/or prevalence values for filtering representative sequences, if desired.

The user can then run the *filtered* mode of the workflow.

```
snakemake --use-conda dada2_pe_denoise --cores 5
# or dada2_se_denoise or deblur_se_denoise

snakemake --use-conda dada2_pe_taxonomy_filtered --cores 5
# or dada2_se_taxonomy_filtered or deblur_se_taxonomy_filtered

snakemake --use-conda dada2_pe_diversity_filtered --cores 5
# or dada2_se_diversity_filtered or deblur_se_diversity_filtered

snakemake --use-conda dada2_pe_report_filtered --cores 5
# or dada2_se_report_filtered or deblur_se_report_filtered
```

Note that the multiple sequence alignment and tree are rebuilt after filtering representative sequences. Depending on the alignment program and parameters, alignment can take several hours to complete.

## Check the output

#### View report

Open your HTML report (`03-reports/report_{method}_{filter}.html`) in [Chrome](https://www.google.com/chrome/) or [Firefox](https://www.mozilla.org/en-US/firefox/new/). To view the linked files: 

* QZV (QIIME 2 visualization): click to download, then drag and drop in [https://view.qiime2.org](https://view.qiime2.org). Empress trees (e.g., `rooted_tree.qzv`) may take more than 10 minutes to load.
* TSV (tab-separated values): click to download, then open in Microsoft Excel or Tabview (command line tool that comes with Tourmaline).
* PDF (portable document format): click to open and view in new tab.

Downloaded files can be deleted after viewing because they are already stored in your Tourmaline directory.

#### Reality check

After completing denoising and core diversity analyses, determine if the results make sense by asking the following questions:

* How many sequences did I start with, and how many are left after denoising?
* Are the representative sequences of similar length or of very different lengths?
* Do the sequence alignment and tree look reasonable?
* Do samples cluster in an expected way in PCoA space?
* Do the taxonomic profiles match expected taxonomic compositions?

## Modify intermediate files

The parameters in `config.yaml` offer many options for customization, but you may want to process your data in a way that Tourmaline doesn't currently support.

You can modify or replace any of the intermediate files in the workflow. As long as the filename and location are the same as the original file, Snakemake will recognize it and use it in the workflow. Here are some examples where this might be useful:

* Alignment: Tourmaline builds the multiple sequence alignment of representative sequences *de novo* (with choice of 3 programs), but you could build it using a different method and then save the alignment as `02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.fasta` and `02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.qza`.
* Tree: Tourmaline builds the tree of representative sequences *de novo* using FastTree, but you could build it using a different method (e.g., with `qiime fragment-insertion`) and then save the tree as `02-output-{method}-{filter}/02-alignment-tree/rooted_tree.qza`.

## Snakemake tips

#### Rules

Snakemake works by executing rules, defined in the `Snakefile`. Rules specify commands and outputs but most critically inputs, which dictate which other rules must be run beforehand to generate those inputs. By defining pseudo-rules at the beginning of the `Snakefile`, we can specify desired endpoint targets as "inputs" that force execution of the whole workflow or just part of it.

When a Snakemake command is run, only those rules that need to be executed to produce the requested target will be run. To make the most of Tourmaline, you might want to familiarize yourself with Snakemake using the [documentation](https://snakemake.readthedocs.io) or follow the tutorial [here](https://github.com/cuttlefishh/tutorials/tree/master/snakemake) to create your own simple Snakemake workflow and understand how it works.

#### Dry run and print shell commands

To see which jobs (rules) and commands will be run by the workflow, use the options `--dryrun` and `--printshellcmds`, respectively. `--dryrun` will prevent the workflow from being executed. `--printshellcmds` can be used with our without `--dryrun`:

```
snakemake --use-conda dada2_pe_denoise --dryrun --printshellcmds
```

#### Regenerate specific files

You can always delete any file you want to regenerate. Then there are several ways to regenerate it:

* Run `snakemake FILE` and Snakemake will determine which rules (commands) need to be run to generate that file.
* Run `snakemake RULE` where the rule generates the desired file as output.

#### Cleanup metadata for file

If a file is being regenerated when you think it shouldn't be, or a symbolic link is not being recognized, it might help to cleanup the metadata of that file. Cleaning up the metadata means that snakemake removes any tracked version info, and any marks that files are incomplete. Do do this, run:

```
snakemake --cleanup-metadata <filenames>
```

#### Directed Acyclic Graph (DAG)

Snakemake provides the command option `--dag` to generate a directed acyclic graph (DAG) of the jobs (rules) that will be run. The DAG is basically a graph that shows the flow and order of rules to reach your desired target, and it can be rendered as a PDF (or other image format) using the program Graphviz on your computer. If you do not have Graphviz installed or are running Tourmaline through Docker, you can install it within the `qiime2-2023.5` conda environment like this:  

```
conda install -c conda-forge graphviz
```

As an example, to get a graph of the rules run when setting the rule `snakemake dada2_pe_report_unfiltered` as your target, run this command:

```
snakemake dada2_pe_report_unfiltered --dag | dot -Tpdf -Grankdir=LR -Gnodesep=0.1 -Granksep=0.1 > dag.pdf
```

For a simplified graph, use the `--rulegraph` option in place of `--dag`:

```
snakemake dada2_pe_report_unfiltered --rulegraph | dot -Tpdf -Grankdir=LR -Gnodesep=0.1 -Granksep=0.1 > rulegraph.pdf
```
