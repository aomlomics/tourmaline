## Configuration

Tourmaline 2 uses three config files, one per step. Example names below reflect defaults; any filename is acceptable. Templates for each file are in the Tourmaline folder.

### 1. QA/QC configuration (Template: config_01_qaqc.yaml)

The qaqc config defines how Tourmaline imports, trims, and summarizes raw reads.

**Core run settings (required)**

```yaml
run_name: my_run                 # unique identifier; used in output folders
output_dir: /path/to/results     # absolute or relative directory for outputs
paired_end: true                 # true for paired-end, false for single-end data
to_trim: true                    # enable primer trimming with Cutadapt
to_merge: false                  # enable vsearch merge (paired-end only)
to_filter: false                # enable feature filtering steps
assay_name: Bacteria-16S-V4V5-Parada  # for metadata reporting
```
Select assay name based on the [NOAA Omics metabarcoding assays](https://github.com/NOAA-Omics/noaa-omics-metabarcoding-assays/blob/main/assays.tsv) controlled vocabulary. If your assay is not available, please create an [issue](https://github.com/NOAA-Omics/noaa-omics-metabarcoding-assays/issues).

**Input data sources (choose one)**

```yaml
raw_fastq_path: /abs/path/to/raw_fastqs       # directory of raw FASTQ(.gz)
trimmed_fastq_path: /abs/path/to/trimmed      # directory of already trimmed FASTQs
sample_manifest_file: 00-data/manifest.tsv    # QIIME 2 manifest (TSV or CSV)
preexisting_fastq_qza: 00-data/demux.qza      # existing demuxed artifact
```

Provide only the fields relevant to your data source; leave others blank. Manifest formats are detailed in [steps/qaqc.md](steps/qaqc.md).

**Primer trimming parameters (required when `to_trim: true`)**

```yaml
fwd_primer: GTGYCAGCMGCCGCGGTAA      # IUPAC supported
rev_primer: GGACTACNVGGGTWTCTAAT
discard_untrimmed: false             # discard reads without primer match
minimum_length: 50                   # post-trimming minimum length (bp)
trimming_threads: 5                  # threads for Cutadapt
```

**Merging + compute options (required when `to_merge: true`)**

```yaml
maxdiffs: 20                         # vsearch merge mismatches
merge_stagger: --p-allowmergestagger # optional vsearch flag
```

### 2. Repseqs configuration (Template: config_02_repseqs.yaml)

Controls ASV generation, filtering, and optional diversity plots.

**Core run settings (required)**

```yaml
run_name: my_run
output_dir: /path/to/results
asv_method: dada2pe          # one of: dada2pe | dada2se | deblur
asv_threads: 5               # threads passed to denoisers
```

**Input data sources (choose one, otherwise will default to** \[output_dir\]/\[my_run-qaqc\])

```yaml
qaqc_run_name: my_qaqc_run         # reuse QA/QC outputs from another run
fastq_qza_file: /abs/path/demux.qza  # external demultiplexed sequences
```

If neither is supplied, the workflow expects demultiplexed reads from the QA/QC step with the same `run_name`.

**Metadata + diversity options**

```yaml
sample_metadata_file: 00-data/metadata.tsv   # optional metadata for summaries/diversity
plot_diversity: true                         # produce alpha/core metrics outputs
alpha_max_depth: 500                         # required when plot_diversity is true
core_sampling_depth: 500                     # required when plot_diversity is true
```

**DADA2 parameters (required when `asv_method` starts with dada2)**

```yaml
dada2_trunc_len_f: 245        # forward truncation length
dada2pe_trunc_len_r: 190      # reverse truncation (paired-end only)
dada2_trim_left_f: 0          # forward trim from left
dada2pe_trim_left_r: 0        # reverse trim from left
dada2_max_ee_f: 2             # forward max expected errors
dada2pe_max_ee_r: 2           # reverse max expected errors
dada2_trunc_q: 2              # truncate at quality score
dada2_pooling_method: pseudo  # independent | pseudo | pooled
dada2_chimera_method: consensus
dada2_min_fold_parent_over_abundance: 1
dada2_n_reads_learn: 1000000
dada2_hashed_feature_ids: --p-hashed-feature-ids  # optional
```

**Deblur parameters (required when `asv_method: deblur`)**

```yaml
deblur_trim_length: 150       # final sequence length (bp)
deblur_trim_left: 0
deblur_mean_error: 0.005
deblur_min_reads: 2
deblur_min_size: 2
deblur_indel_max: 3
reference_seqs: 00-data/ref.qza  # required reference set
```

**Post-denoising filtering (required if `to_filter` is `True`)**

```yaml
repseq_min_length: 0
repseq_max_length: 0
repseq_min_abundance: 0
repseq_min_prevalence: 0
repseq_min_frequency: 0
repseq_min_samples: 0
```

### 3. Taxonomy configuration (Template: config_03_taxonomy.yaml)

Defines how representative sequences are assigned taxonomy and summarized.

**Core run settings (required)**

```yaml
run_name: my_run
output_dir: /path/to/results
classify_method: naive-bayes      # options: naive-bayes | consensus-blast | consensus-vsearch | bt2-blca | revamp
taxa_ranks: kingdom,phylum,class,order,family,genus,species
collapse_taxalevel: 7             # taxonomy level for collapsed table
classify_threads: 10
```

**Input data sources (choose one)**

```yaml
repseqs_run_name: my_repseqs_run          # reuse outputs from another run
repseqs_qza_file: /abs/path/repseqs.qza   # external representative sequences
table_qza_file: /abs/path/table.qza       # external feature table
```

If no external inputs are supplied, the workflow uses artifacts produced by the Repseqs step with matching `run_name`.

**Reference database parameters**

```yaml
database_name: silva-138_1
refseqs_file: 00-data/silva-seqs.qza    # required unless using pretrained classifier
taxa_file: 00-data/silva-tax.qza        # required unless using pretrained classifier
sample_metadata_file: 00-data/metadata.tsv  # optional for barplots
```

**Naive Bayes options**

```yaml
pretrained_classifier: /abs/path/classifier.qza  # optional, overrides refseqs/taxa files
skl_confidence: 0.7                              # confidence threshold
```

**Consensus BLAST/VSEARCH options**

```yaml
perc_identity: 0.8
query_cov: 0.8
min_consensus: 0.51
```

**BT2-BLCA options**

```yaml
bowtie_database: /abs/path/bowtie2_index/   # optional; auto-built if omitted
confidence_thres: 0.8
```

**REVAMP options**

```yaml
revamp_dir: /abs/path/REVAMP                        # REVAMP clone
revamp_blastdb: /abs/path/blastdb                   # nt volumes plus prepared taxdump/
revamp_blast_results: /abs/path/ASV_blastn_nt.btab  # optional; BLAST run elsewhere
revamp_blast_mode: mostEnvOUT                       # allIN | allEnvOUT | mostEnvOUT
revamp_query_cov: 90                                # percent of ASV length a hit must cover
revamp_taxonomy_cutoffs: "97,95,90,80,70,60"        # percent ID cutoffs, ordered S,G,F,O,C,P
```

Used only when `classify_method` is `revamp`, which needs the `revamp` conda environment
and ignores `refseqs_file` / `taxa_file` / `pretrained_classifier`. `revamp_blast_mode`
applies when Tourmaline runs BLAST itself; with `revamp_blast_results` supplied it is
recorded but not applied. Suggested cutoffs are `97,95,90,80,70,60` for rRNA genes and
`95,92,87,77,67,60` for protein-coding genes. See
[Taxonomy step](steps/taxonomy.md#revamp) for database preparation, running BLAST on
another machine, and how REVAMP's output differs from the other methods.

**Krona plot options (any classify method)**

```yaml
make_krona: False        # build figures/{run_name}-krona.html
krona_per_sample: True   # add one Krona dataset per sample
```

Off by default because it needs a `krona` conda environment
(`conda create -c conda-forge -c bioconda -n krona krona`). Configs written before this
option was added still work; the plot is simply not built. See
[Taxonomy step](steps/taxonomy.md#krona-plots).

**Additional classifier flags**

```yaml
classify_params: --verbose   # appended to the chosen classifier command
```

See [Running](running.md) for multi-step invocation and [External Data](external_data.md) for conversions and artifact preparation tips.

### 4. Tax-credit configuration (Template: config_04_tax_credit.yaml)

The tax-credit step benchmarks reference databases using [tax-credit](../tax-credit/) simulations and Tourmaline taxonomy assignment rules. It does not require outputs from the repseqs or taxonomy steps.

**Core run settings**

```yaml
run_name: mifish_tax_credit
output_dir: ../v2-results
tax_credit_package_dir: ../tax-credit   # pip install -e this path in qiime2 env
```

**Reference databases (one or more)**

```yaml
reference_databases:
  - id: my_database
    refseqs_file: /path/to/sequences.fasta   # or .qza
    taxa_file: /path/to/taxonomy.txt         # or .qza
    fwd_primer: GCCGGTAAAACTCGTGCCAGC
    rev_primer: CATAGTGGGGTATCTAATCCCAGTTTG
    fwd_primer_id: MiFishF
    rev_primer_id: MiFishR
    read_length: 250            # required when truncate is true
    min_read_length: 140
    trim_primers: true          # false for pre-trimmed amplicon references (e.g. rCRUX)
    truncate: true
    pretrained_classifier:      # optional, mock-community naive-bayes only
```

Primers and the simulation settings are only needed by the simulated evaluation methods.

**Evaluation methods (one or more)**

```yaml
evaluation_methods:
  - cross-validated          # taxonomy-aware CV folds
  - cross-validated-trad     # traditional KFold CV
  - novel-taxa               # novel-taxa simulation
  - self-validated           # full database classified against itself
  # - mock-community         # see Mock community evaluation below
```

**Simulation parameters**

```yaml
iterations: 10
novel_taxa_levels: [6, 5, 4, 3]
cv_recall_max_level: 6
novel_recall_min_level: 3
force_regenerate: false
```

In novel-taxa folds the expected taxonomy of each query is truncated to the deepest rank still present in that fold's reference. Removing the novel taxon can remove its parent as well (a monotypic genus goes with its only species), and a rank the reference no longer holds cannot be returned by any classifier, so grading against it would count every method as wrong whatever it did. Queries with no rank left in the reference are dropped from the fold. Dataset generation reports how many queries this affects per fold; novel-taxa scores are therefore not comparable with runs generated before this behaviour.

`cv_recall_max_level` controls cross-validated assignment manifest generation (default `6`; `min_level` is always `max_level - 1` so each CV fold is listed once). Per-database simulation settings (`read_length`, `min_read_length`, `trim_primers`, `truncate`) are configured on each `reference_databases` entry (see above).

**Taxonomic assignment** — uses the same keys as the taxonomy step (`classify_method`, `skl_confidence`, `classify_params`, etc.). Assignment runs via Snakemake rules shared with `taxonomy_step.Snakefile`, not tax-credit shell templates.

```yaml
classify_method: naive-bayes
classify_threads: 5
nb_confidence_values: [0.7]
blca_confidence_values: [0.8]
blca_perc_identity: 0.8
blca_query_cov: 0.8
skl_confidence: 0.7
confidence_thres: 0.8
fit_params: "--p-feat-ext--ngram-range '[7,7]' --p-classify--alpha 0.001"
generate_plots: true
```

`classify_methods` accepts `naive-bayes`, `consensus-blast`, `consensus-vsearch`, `bt2-blca` and `revamp`. Assignment jobs write method-relevant parameters to `assignment_manifest.tsv`; unused fields are left blank. List-valued `perc_identity`, `query_cov`, and `min_consensus` expand into parameter sweeps for consensus methods. bt2-blca has its own cutoffs, `blca_perc_identity` (BLCA `-b`) and `blca_query_cov` (BLCA `-l`, minimum hit length relative to the query), which also accept lists and do not read the consensus keys. Both are required whenever `bt2-blca` is in `classify_methods`: a config without them (such as an older config that relied on `perc_identity` / `query_cov` for bt2-blca) stops with an error before any work starts.

**REVAMP options (mock-community only)**

```yaml
revamp_dir: /abs/path/REVAMP            # REVAMP clone
revamp_blastdb: /abs/path/blastdb       # NCBI nt directory with a prepared taxdump/
revamp_query_cov_values: [90]           # percent of ASV length a BLAST hit must cover
revamp_taxonomy_cutoffs_values:         # percent identity cutoffs, ordered S,G,F,O,C,P
  - "97,95,90,80,70,60"
  - "95,92,87,77,67,60"
revamp_blast_mode_values: [mostEnvOUT]  # allIN | allEnvOUT | mostEnvOUT
```

`revamp` classifies against a local NCBI `nt` BLAST database, which is far too large to simulate cross-validated, novel-taxa or self-validated datasets from, so it runs **only for `mock-community`**. Listing it alongside other evaluation methods is fine — it is skipped for them with a message, and the other classify methods still run. Listing it with no `mock-community` in `evaluation_methods` is an error.

Because `nt` is not a QIIME reference artifact, it is declared as a `reference_databases` entry carrying `revamp: true` and no `refseqs_file` / `taxa_file`:

```yaml
reference_databases:
  - id: ncbi-nt
    revamp: true
```

Only `revamp` classifies against that entry, and `revamp` classifies against nothing else — every other database/method pairing is skipped. Its expected set must use the NCBI taxonomy backbone, since that is what REVAMP assigns; the backbone check is skipped for it because there is no local reference taxonomy to compare against.

All three sweep keys accept a single value or a list. `revamp_query_cov_values` and `revamp_taxonomy_cutoffs_values` are applied to BLAST output, so their combinations reuse one BLAST search — the same fit-job sharing bt2-blca uses for its bowtie2 index. `revamp_blast_mode_values` changes the search itself, so each mode costs another full BLAST against `nt`; with precomputed results it is ignored (the filtering already happened) and a message says so. Parameter sets appear in summaries and plots as, for example, `qc90-cut97_95_90_80_70_60-mostEnvOUT`.

BLASTing `nt` usually happens on the machine that holds it. Give each mock dataset a `blast_results` btab and no BLAST runs locally; leave it empty and one BLAST job runs per dataset and mode, shared by every parameter combination. See [Taxonomy step](steps/taxonomy.md#revamp) for the BLAST command, the `revamp` conda environment and database preparation.

**Plotting and log analysis**

```yaml
generate_plots: true
plot_types: [boxplot, pointplot, heatmap, stacked_bar, best_run_stacked_bar]
plot_metrics: [Precision, Recall, F-measure, match_ratio, underclassification_ratio, overclassification_ratio, misclassification_ratio]
plot_ranks: [genus, species]
best_run_rank: species
generate_log_analysis: true
log_analysis_ranks: [species, genus, family]
```

Cross-validated and self-validated metrics are plotted per fold and per taxonomic rank. Precision, Recall and F-measure come from the evaluation summary; the four classification ratios are recomputed from each fold's `classification_accuracy_log.tsv`. Pointplots show each metric across ranks, boxplots show the spread with one panel per rank in `plot_ranks` (optional, default `[genus, species]`). Both group by dataset and method, pooling all folds and parameter sets of a method; heatmaps show parameter sets separately. CSVs keep numeric levels (1 = phylum … 6 = species); plot axes use rank names.

Plot files are named `<evaluation method>-<metric>-<plot type>.pdf` inside each evaluation method's folder, e.g. `plots/cross-validated/cross-validated-F-measure-heatmap.pdf`, and are read from the summary named in `summary_filenames`. Each classify method keeps the same colourblind-safe colour in every plot unless `plot_color_palette` sets a seaborn palette name or a method → colour mapping. Axes and heatmap colour scales zoom (and say so) when every value is within 0.25 of 0 or 1, so near-zero novel-taxa scores stay visible. `<evaluation method>-classification-ratios-stacked-barplot.pdf` has one row per method + parameter set and one column per reference database. PDF text is embedded as editable TrueType (Arial, or DejaVu Sans where Arial is not installed).

`best_run_stacked_bar` picks, for each reference database and each metric in `plot_metrics`, the method + parameter combination with the best score averaged over folds (lowest for the mis-, over- and underclassification ratios, highest otherwise; ties go to the first run alphabetically). Cross-validated and self-validated runs are compared at `best_run_rank` (optional, default `species`); novel-taxa runs are compared within each novel level. It writes one figure per evaluation method, `<evaluation method>-best-run-stacked-barplot.pdf`, holding every database in that evaluation: databases (and novel levels) are rows, metrics are columns, and each panel is a stacked barplot of classification ratios by rank for that metric's winning run. It lists the selections (score, fold count, number of tied runs) in `summaries/best_runs/<evaluation method>/<summary>.csv`.

Novel-taxa results are split by novel level (`L5`, `L6`, …). Stacked barplots are written per level (`novel-taxa-classification-ratios-stacked-barplot-L6.pdf`) and show only ranks above the novel rank. Method/parameter sensitivity is computed at the rank just above the novel rank (`novel-taxa-method-parameter-misclassified-L6-genus.pdf`) instead of `log_analysis_ranks`.

Log-analysis tables (`taxon_error_profiles.csv`, `confusion_pairs.csv`, `cross_fold_stability.csv`, `method_parameter_misclassified-*.csv`) include an `expected_rank` column: the deepest rank named in the expected taxonomy. It is shallower than the analysis rank when the reference lacks that rank or, in cross-validated folds, when the taxon's lower ranks are absent from the training fold. Misclassification heatmaps (`<evaluation method>-method-parameter-misclassified-<rank>.pdf`) show only taxa whose `expected_rank` is the plotted rank, with one panel per reference database holding up to `log_analysis_top_n` taxa ranked by their highest misclassification in any run (taxa never misclassified are left out); hatched cells mean the taxon had no reads in that run.

#### Mock community evaluation

Add `mock-community` to `evaluation_methods` to score taxonomy assignments of real sequencing data from mock communities, whose make-up you know. Every classify method and parameter set in the config is run on the mock ASVs against every reference database assigned to an expected set, and each result is compared with what the mock should contain. Mock evaluation needs no simulation, so `fwd_primer` / `rev_primer` are optional for databases used only here.

**Inputs**

| Input | Config key | Required | Format |
|---|---|---|---|
| Read counts | `datasets[].feature_table` | yes | TSV with ASVs as rows and samples as columns (first column = feature id), BIOM, or a repseqs `table.qza` |
| ASV sequences | `datasets[].rep_seqs` | yes | FASTA or a repseqs `repseqs.qza`; ids must match the feature table |
| Expected composition | `expected_sets[].composition` | this, `asv_taxonomy`, or both | TSV, see below |
| Known ASV taxonomy | `expected_sets[].asv_taxonomy` | this, `composition`, or both | TSV, see below |

Expected composition: one row per unique taxonomy, one column per mock sample, relative abundance as values. The first column header can be anything. Missing values are 0, and each sample column is rescaled to sum to 1. Duplicate rows (after the cleanup described under *Taxonomy strings*) are summed.

```
Taxonomy	mock-even	mock-staggered
Eukaryota;Chordata;Actinopteri;Lophiiformes;Lophiidae;Lophius;Lophius americanus	0.25	0.60
Eukaryota;Chordata;Actinopteri;Lampriformes;Lampridae;Lampris;Lampris guttatus	0.25	0.30
Eukaryota;Chordata;Chondrichthyes;Hexanchiformes;Hexanchidae;Hexanchus	0.50	0.10
```

Known ASV taxonomy: the true taxonomy of each ASV, with headers `Feature ID` and `Taxon` (the QIIME 2 `taxonomy.tsv` layout; extra columns are ignored). ASVs left out, or with an empty or `Unassigned` taxon, count as unknown.

```
Feature ID	Taxon
672935cdade199c3f4433afc5ba327fd	Eukaryota;Chordata;Actinopteri;Lophiiformes;Lophiidae;Lophius;Lophius americanus
```

What each expected set enables:

| Provided | Mock samples | Metrics |
|---|---|---|
| `composition` only | composition columns that are in the feature table | TAR, TDR, Bray-Curtis |
| `asv_taxonomy` only | feature-table samples with reads from ASVs of known taxonomy; the expected composition is built from those reads | TAR, TDR, Bray-Curtis, precision, recall, F-measure |
| both | composition columns that are in the feature table | TAR, TDR and Bray-Curtis use `composition`; precision, recall and F-measure use `asv_taxonomy` |

Feature-table samples that are not mock samples (blanks, field samples) are ignored and listed, with the reason, in `data/mock-community/datasets/<dataset>/excluded_samples.tsv`, so you can point `feature_table` at a whole run. Set `datasets[].samples` to evaluate only some samples.

**Taxonomic backbones: one expected set per backbone**

Reference databases can name and place the same organism differently. For example, *Antigonia combatia* is in `Acanthuriformes;Antigoniidae` in one MiFish database and `Caproiformes;Caproidae` in another. Expected taxa written against one backbone would be scored as errors against a database that uses the other. Each expected set therefore lists the `databases` it applies to, and must use those databases' taxonomy. A database can belong to only one expected set. Databases in no expected set are not used for mock evaluation (they can still be used by the other evaluation methods).

Before any assignment runs, the datasets phase checks every expected taxonomy against each database in its set, at each rank in `eval_ranks`, and writes `data/mock-community/expected/<set>/backbone_check-<database>.tsv`:

| status | meaning |
|---|---|
| `found` | the taxonomy, down to that rank, exists in the database |
| `different_lineage` | the name at that rank exists in the database under a different lineage (listed in `database_lineages`); almost always a backbone mismatch |
| `not_in_database` | the name is not in the database at that rank; the database may lack the taxon, or it is named differently |

It also warns when the number of ranks or the rank-prefix style (`g__Name` vs `Name`) differs between the expected taxa and the database. With `backbone_check: warn` (default) the run continues and a per-rank count is printed. With `backbone_check: error` any taxon that isn't `found` stops the run.

**Changing datasets or databases after a first run.** The evaluation plan is written by the staging phase, which does not re-run when only the config changes. If you add a mock dataset or a reference database (for example an `ncbi-nt` entry for revamp) to an existing run name, delete the `.datasets.done` marker in the run output directory — or pass `--forcerun tax_credit_prepare_datasets` — so staging rebuilds the plan. The manifest phase compares the config against `data/mock-community/staged_inputs.tsv` and stops with this message if they disagree, rather than quietly skipping the new database.

**Config**

```yaml
mock_community:
  datasets:                       # one or more sequencing datasets with mock samples
    - id: pmel_mifish_mocks
      feature_table: /path/to/table.tsv     # .tsv, .biom or table.qza
      rep_seqs: /path/to/asvs.fasta          # .fasta or repseqs.qza
      samples:                               # optional: only these samples
      blast_results: /path/to/mock.btab      # optional, revamp only: BLASTn vs nt
  expected_sets:                  # one per taxonomic backbone
    - id: gomex_backbone
      databases: [mifishGom]                 # reference_databases ids using this backbone
      composition: /path/to/composition_gomex.tsv
      asv_taxonomy:                          # optional; enables precision/recall/F-measure
      ranks:                                 # optional; default taxa_ranks
    - id: addJ_backbone
      databases: [addJ]
      composition:
      asv_taxonomy: /path/to/asv_taxonomy_addJ.tsv
      ranks:
  eval_ranks: [family, genus, species]      # ranks to score (default family, genus, species)
  min_relative_abundance: 0                  # observed taxa must exceed this fraction of reads for TAR/TDR
  legacy_unresolved_taxa: false              # true: original tax-credit TAR/TDR (see below)
  backbone_check: warn                       # warn | error
  plot_metrics:                              # default: TAR, TDR, Bray-Curtis, Precision, Recall, F-measure
  best_run_metric: auto                      # metric that picks each method's best parameter set
  composition_top_n: 12                      # taxa coloured individually in composition plots
```

`ranks` names the positions in the expected set's taxonomy strings (default: the top-level `taxa_ranks`), so `eval_ranks` can be given by name. A reference database entry can also set `pretrained_classifier: /path/to/classifier.qza`. naive-bayes then classifies the mock ASVs with that classifier instead of fitting one on the full database; the parameter id is `nb-pretrained-conf<confidence>`.

**Metrics** (per mock sample, per rank, per reference database, method and parameter set)

At each rank, taxonomy strings are cut after that rank. A string is *resolved* at a rank if it has a name there, so an ASV assigned only to genus, and an expected taxon only known to genus, are unresolved at species.

- **Taxon Accuracy Rate (TAR)**: the fraction of observed taxa that are expected, |observed ∩ expected| / |observed|. Observed taxa have more than `min_relative_abundance` of the sample's reads; expected taxa have expected abundance above 0. Only resolved taxa count on either side. NaN when no taxa are observed.
- **Taxon Detection Rate (TDR)**: the fraction of expected taxa that are observed, |observed ∩ expected| / |expected|. NaN when no taxa are expected at that rank.
- **Bray-Curtis**: dissimilarity between the expected composition and the observed read proportions at that rank (0 = identical, 1 = nothing shared; lower is better). Unresolved assignments keep their truncated taxonomy, and unassigned reads are pooled as `Unassigned`, so reads that are not assigned still count. No abundance threshold is applied. PCR and copy-number bias mean observed proportions rarely match the true mix even with perfect classification, so use Bray-Curtis to compare methods rather than as an absolute score.
- **Precision, Recall, F-measure** (needs `asv_taxonomy`): each ASV's assignment is compared with its known taxonomy and weighted by its reads in the sample. A match is a true positive; underclassification (a correct but shallower assignment) is a false negative; overclassification and misclassification are a false positive and a false negative. This is the scoring used for the simulated evaluation methods. `ASV Precision`, `ASV Recall` and `ASV F-measure` weight every ASV equally. `match_ratio`, `underclassification_ratio`, `overclassification_ratio` and `misclassification_ratio` are the read-weighted fraction of each outcome. ASVs of unknown taxonomy are left out; `reads_scored_fraction` shows how many reads were scored.

Metrics that cannot be computed are left empty (NaN), not `-1`.

**Differences from the original tax-credit mock evaluation**

- *Unresolved taxa in TAR/TDR.* Original tax-credit treats a truncated string such as `…;Lophiidae;Lophius` as its own taxon at species rank, so an assignment that stops at genus lowers TAR as a false taxon, and an expected taxon known only to genus is counted as expected at species. Here both are left out of TAR/TDR at ranks where they are unresolved (a shallow assignment is not a claim of a wrong species), and they are scored by Bray-Curtis and recall instead. Set `legacy_unresolved_taxa: true` to count them as taxa, including `Unassigned`, as in the original. Scores are only comparable with published tax-credit results in legacy mode.
- *No-detection cases.* The original reports TAR and TDR as 0 when no expected taxon is observed, even if nothing was observed. Here TAR is NaN when nothing is observed.
- *Inputs.* Expected results come from one composition table and/or one ASV taxonomy table per backbone, not per-dataset `expected/` directories with BIOM tables collapsed per level, `expected-taxonomy.tsv` and `trueish-taxonomies.tsv`. Ranks are named, not numbered.
- *Precision and recall.* Weighted by reads as in the original (`per_seq_precision`), but computed whenever `asv_taxonomy` is given, alongside unweighted per-ASV scores and classification ratios.
- *Recomputation.* Snakemake scores each assignment job separately and rescores only jobs whose assignments changed; the `force` / `append` / `backup` options are gone.
- *Taxonomy strings.* Whitespace around ranks is removed, trailing `NA`/empty ranks are dropped (internal `NA` ranks keep their place), and `Unassigned`, `Unclassified` and `No blast hit` all mean unassigned. The original's removal of `[]()` characters is not done.

**Outputs**

```
[run_name]-tax-credit/
├── data/
│   ├── ref_dbs/<database>/ref_seqs.qza, ref_taxa.qza
│   ├── mock-community/
│   │   ├── dataset_log.txt                   # staging messages: kept/excluded samples, unscored ASVs, backbone check
│   │   ├── evaluation_plan.tsv               # dataset × expected set × database × samples
│   │   ├── datasets/<dataset>/
│   │   │   ├── rep_seqs.qza, feature_table.tsv
│   │   │   └── excluded_samples.tsv          # samples not evaluated, per expected set, with reason
│   │   └── expected/<set>/
│   │       ├── composition.tsv, asv_taxonomy.tsv
│   │       ├── backbone_check-<database>.tsv
│   │       └── unscored_asvs-<dataset>.tsv   # ASVs left out of precision/recall (asv_taxonomy sets only)
│   └── results/mock-community/<dataset>/<database>/<method>/<parameters>/taxonomy.tsv
├── summaries/
│   ├── mock_community_metrics.tsv            # one row per sample × rank × run
│   ├── mock_community_composition.tsv        # expected vs observed abundance per taxon
│   ├── mock-community/per-job/                # per-run scores the summaries are built from
│   └── best_runs/mock-community/
└── plots/mock-community/
```

`dataset_log.txt` is rewritten each time inputs are staged (the datasets phase) and records everything staging reports: feature and sample counts per dataset, the samples evaluated and excluded for each expected set, features with no sequence in `rep_seqs`, unscored ASVs, the backbone check counts and warnings, and the error that stopped staging, if any. The same lines are printed to the terminal with a `[mock-community]` prefix. Messages from later scoring jobs (features with no assignment in a run) are only printed to the terminal.

`excluded_samples.tsv` has one row per excluded sample and expected set (`sample_id`, `expected_set`, `reason`). Reasons: `not in composition`; `in composition but not in feature table`; `no reads from ASVs with known taxonomy` (sets with only `asv_taxonomy`); `not in datasets.samples`. The file has only a header when nothing was excluded.

`unscored_asvs-<dataset>.tsv` lists feature-table ASVs without a known taxonomy (`Feature ID`, `reason` = `not in asv_taxonomy` or `empty or Unassigned in asv_taxonomy`, `reads_in_mock_samples`, `reads_total`), sorted by reads in the evaluated mock samples. It is written only for expected sets with `asv_taxonomy` and pairs that have mock samples.

`mock_community_metrics.tsv` columns: `MockDataset`, `Reference`, `ExpectedSet`, `Method`, `Parameters`, `SampleID`, `rank`, `level`, the metrics above, and the counts `n_expected_taxa`, `n_observed_taxa`, `n_shared_taxa`, `n_asvs_scored`, `reads_scored_fraction`. `mock_community_composition.tsv` has one row per run, sample, rank and taxon, with `expected` and `observed` relative abundance and whether the taxon is `resolved` at that rank.

Plots are drawn per mock dataset, for the metrics in `mock_community.plot_metrics` and the `plot_types` in the config:

- `boxplot`: one point per sample and parameter set, by reference database and method, with one panel per rank in `plot_ranks` that is also in `eval_ranks`. A second figure, `…-best-boxplot.pdf`, shows only each method's best parameter set, so one point per sample.
- `pointplot`: metric by rank, with `…-best-pointplot.pdf` again restricted to the best parameter sets.
- `heatmap`: mean over samples for each method + parameter set, by database and rank.
- `stacked_bar`: read-weighted classification ratios by rank (needs `asv_taxonomy`).
- `best_run_stacked_bar`: `mock-community-<dataset>-composition-<rank>.pdf` draws the expected composition next to the observed composition of each method's best parameter set. One figure per mock dataset holds every reference database: databases are rows, mock samples are columns, and all panels share the taxon colours. `summaries/best_runs/mock-community/mock_community_best_runs.csv` also lists the best run per database for every metric in `plot_metrics`.

**Choosing the best run.** `mock_community.best_run_metric` sets the metric that picks each method's best parameter set for the composition plot and the best-run box and point plots. Runs are compared at `best_run_rank` (or the deepest `eval_ranks` entry when `best_run_rank` is not one of them), averaged over the mock samples; Bray-Curtis is minimised and every other metric maximised. The default `auto` uses F-measure for databases whose expected set has an `asv_taxonomy` and Bray-Curtis for the rest, so a run holding both kinds of expected set picks per database. Naming a metric (for example `Taxon Detection Rate`) uses it everywhere; a database with no score for that metric gets no best run and a message says so. The picks are written to `summaries/best_runs/mock-community/mock_community_best_per_method.csv` with the metric, direction, value and number of tied runs.

Log analysis (`generate_log_analysis`) applies to the simulated evaluation methods only.

**Smoke test**: `config_04_tax_credit_test_mock.yaml` runs every mock code path on the small test databases and the fixture in `00-data/tax-credit-test/mock/` (built by `scripts/make_mock_test_fixture.py`). The two test databases disagree on the backbone for a few species, so it uses one expected set for each.

Install tax-credit in the QIIME 2 amplicon environment before running:

```bash
pip install -e ../tax-credit
```


