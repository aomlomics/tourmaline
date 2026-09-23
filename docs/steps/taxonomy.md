## Taxonomy Step

Assigns taxonomy and generates visualizations using one of five methods:

- Naive Bayes (QIIME 2 classify-sklearn)
- Consensus BLAST
- Consensus VSEARCH
- BT2-BLCA (Anacapa)
- REVAMP (BLASTn against a local NCBI nt database; see [REVAMP](#revamp) below)

### Configuration

```yaml
classify_method: naive-bayes        # or consensus-blast, consensus-vsearch, bt2-blca, revamp
classify_threads: 4
collapse_taxalevel: 0               # build collapsed table at this rank (0 disables)
taxa_ranks: kingdom,phylum,class,order,family,genus,species

# Inputs from Repseqs step or external artifacts (choose one source)
repseqs_run_name: repseqs_run        # reuse artifacts from another repseqs run
repseqs_qza_file: /abs/path/repseqs.qza   # externally supplied representative sequences
table_qza_file: /abs/path/table.qza       # externally supplied feature table

# Reference database (required unless using a pretrained classifier)
database_name: PR2
refseqs_file: 00-data/refseqs.fna         # FeatureData[Sequence] artifact or FASTA
taxa_file: 00-data/reftax.tsv             # FeatureData[Taxonomy] artifact or TSV
pretrained_classifier: /abs/path/classifier.qza  # optional for naive-bayes
bowtie_database: /abs/path/bt2/index/            # optional cache for bt2-blca

# Thresholds / method options
skl_confidence: 0.7          # naive-bayes confidence cutoff
perc_identity: 0.8           # consensus BLAST/VSEARCH identity threshold
query_cov: 0.8               # consensus BLAST/VSEARCH query coverage
min_consensus: 0.51          # consensus BLAST/VSEARCH agreement fraction
confidence_thres: 0.8        # bt2-blca confidence cutoff
classify_params: --verbose   # optional extra args for the classifier

# REVAMP options (classify_method: revamp only)
revamp_dir: /abs/path/REVAMP              # REVAMP clone
revamp_blastdb: /abs/path/blastdb         # nt volumes plus a prepared taxdump/
revamp_blast_results: /abs/path/ASV_blastn_nt.btab  # optional: BLAST run elsewhere
revamp_blast_mode: mostEnvOUT             # allIN | allEnvOUT | mostEnvOUT
revamp_query_cov: 90                      # percent of ASV length a hit must cover
revamp_taxonomy_cutoffs: "97,95,90,80,70,60"  # percent ID cutoffs, ordered S,G,F,O,C,P

# Krona plot (optional, works with any classify method)
make_krona: False            # requires the `krona` conda env
krona_per_sample: True       # one dataset per sample alongside the all-samples plot

# Optional metadata for barplots
sample_metadata_file: 00-data/metadata.tsv
```

### Outputs

- Taxonomy assignments and visualizations (barplots, summaries)
- Optional collapsed table at `collapse_taxalevel`
- Combined ASV/taxonomy/sequence TSV for downstream FAIR metadata workflows
- Exported taxonomy TSV alongside the QIIME 2 taxonomy artifact
- Optional Krona plot at `figures/{run_name}-krona.html` when `make_krona` is `True`

### Krona plots

Setting `make_krona: True` adds an interactive [Krona](https://github.com/marbl/Krona/wiki)
plot for any classify method. `scripts/taxonomy_to_krona.py` turns the taxonomy artifact
and feature table into Krona text files under `figures/krona_inputs/`, and `ktImportText`
renders them into one self-contained HTML file. Rank prefixes (`k__`, `d__`) are stripped,
trailing empty or `NA` ranks are dropped, and unassigned features are grouped under
`Unassigned`.

The plot always includes an `all_samples` dataset summed across the run; with
`krona_per_sample: True` each sample is added as its own dataset, selectable from the
dropdown at the top of the plot. Set it to `False` for runs with many samples.

Requires a `krona` environment:

```bash
conda create -c conda-forge -c bioconda -n krona krona
```

No Krona taxonomy database is needed — `ktUpdateTaxonomy.sh` applies only to Krona's
taxid-based importers, and Tourmaline supplies resolved lineage strings instead. The
result is a plain HTML file rather than a QIIME 2 visualization, so open it directly in a
browser rather than through `qiime tools view`.

### REVAMP

[REVAMP](https://github.com/McAllister-NOAA/REVAMP) assigns taxonomy by BLASTing ASVs
against a local NCBI `nt` database and merging all best hits to their lowest common
ancestor, limited by percent identity cutoffs per rank. Tourmaline calls REVAMP's
taxonomy scripts directly through `scripts/run_revamp_taxonomy.sh` rather than running
`revamp.sh`, which also performs read trimming, denoising, tables and figures. Nothing
inside the REVAMP clone is modified, so it can be updated with `git pull`.

`refseqs_file`, `taxa_file` and `pretrained_classifier` are unused: `nt` is the reference
database and NCBI `taxonomy` the reference taxonomy. `taxa_ranks` must list seven ranks,
because REVAMP always assigns kingdom through species.

**Requirements**

- A `revamp` conda environment:
  ```bash
  mamba create -n revamp -c conda-forge -c bioconda \
    "blast>=2.13" "taxonkit>=0.20" r-base r-dplyr bioconductor-biostrings \
    perl perl-list-moreutils krona
  ```
  `taxonkit` must be 0.20 or newer for its `reformat2` command. Krona is optional; without
  it REVAMP skips its KRONA HTML plots and assignments are unaffected.
- A database directory (`revamp_blastdb`) containing the `nt` volumes and a `taxdump/`
  subdirectory prepared by REVAMP's `ncbi_db_cleanup.sh`: `names.dmp`, `nodes.dmp`,
  `merged.dmp`, `delnodes.dmp`, `common_names.dmp`, and the taxid exclusion lists used by
  `allEnvOUT` / `mostEnvOUT`. `TAXONKIT_DB` is pointed at this directory for the duration
  of the run, so no files are copied into the user's own taxonkit directory.

**Running BLAST separately.** The `nt` database often lives on another machine. Run BLAST
there on the ASVs Tourmaline exported to
`{run_name}-taxonomy/revamp/dada2/ASVs.fa`, then point `revamp_blast_results` at the
result. The column order matters, since REVAMP reads the fields positionally:

```bash
blastn -db <blastdb>/nt -query ASVs.fa \
  -outfmt '6 qseqid pident length staxids sacc' \
  -subject_besthit -max_target_seqs 4000 -num_threads <N> \
  -negative_taxidlist <blastdb>/taxdump/taxid_exclusion_list_leavesinUnclassified.txt \
  -out ASV_blastn_nt.btab
```

Taxid filtering additionally needs BLAST's `taxdb` files (`taxdb.btd`, `taxdb.bti`,
`taxonomy4blast.sqlite3`) on the `BLASTDB` search path; without them BLAST reports that
the option "requires additional data files". Drop `-negative_taxidlist` and set
`revamp_blast_mode: allIN` if they are unavailable.

The same applies when Tourmaline runs BLAST itself: `allEnvOUT` and `mostEnvOUT` filter by
taxid, so `run_revamp_taxonomy.sh` checks for the `taxdb` files before starting and stops
with installation instructions if they are missing. This check matters because BLAST
otherwise prints that message and **keeps going with the exclusion list unapplied**,
producing `allIN` results labelled as `mostEnvOUT` after a full-length search. Install
them into the database directory with:

```bash
cd /path/to/blastdb && update_blastdb.pl taxdb && tar -xzf taxdb.tar.gz
```

The database directory is prepended to `BLASTDB` for the search, so the files are found
whether they sit beside the `nt` volumes or elsewhere on that path.

With supplied results, only `taxdump/` is needed locally — the `nt` volumes are not read.
The run fails if the feature table and representative sequences describe different ASVs,
or if the BLAST file contains query IDs absent from the representative sequences. ASVs
with no BLAST hit are reported as a count and come back `Unassigned`.

**Output differences.** The taxonomy TSV's third column is `percent_id`: the best-hit
BLAST percent identity on a 0-100 scale, since REVAMP produces no confidence score. Other
methods write `Confidence` on a 0-1 scale, and this column flows through to
`-asv_taxa_features.tsv`. Taxon strings are plain NCBI names with no rank prefixes. Ranks
REVAMP could not resolve are `NA`; trailing `NA` ranks are dropped and internal ones keep
their position, so a genus never shifts into the family slot. REVAMP's own intermediate
files — formatted BLAST hits, taxonkit lineages, common names, KRONA inputs, per-ASV
taxonomy table — are kept under `{run_name}-taxonomy/revamp/`.

See [External Data](../external_data.md) for conversions if starting outside Tourmaline.

