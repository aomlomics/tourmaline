## Taxonomy Step

Assigns taxonomy and generates visualizations using one of five methods:

- Naive Bayes (QIIME 2 classify-sklearn)
- Consensus BLAST
- Consensus VSEARCH
- BT2-BLCA (Anacapa)
- REVAMP (BLASTn against a local NCBI nt database; see [REVAMP](#revamp) below)
- REVAMP with your own reference database and taxonomy, any number of ranks
  (`revamp-custom`; see [below](#revamp-with-a-custom-reference-database-revamp-custom))

### Configuration

```yaml
classify_method: naive-bayes        # or consensus-blast, consensus-vsearch, bt2-blca, revamp, revamp-custom
classify_threads: 4
collapse_taxalevel: 7               # rank (1 to number of taxa_ranks) for the additional collapsed count table
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
max_accepts: 10              # consensus BLAST/VSEARCH max hits per query (`all` for vsearch)
blast_database:              # consensus-blast: prebuilt BLAST db dir; enables classify_threads
build_blast_database: false  # consensus-blast: build one from refseqs_file instead
confidence_thres: 0.8        # bt2-blca confidence cutoff
blca_bowtie2_max_hits: 100   # bt2-blca: bowtie2 -k, alignments kept per query
blca_bootstraps: 100         # bt2-blca: BLCA -n, bootstrap replicates per query
blca_muscle_maxiters: 16     # bt2-blca: muscle -maxiters per query alignment
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

### Reference file formats

`refseqs_file` must be a QIIME 2 artifact (`.qza`) or FASTA (`.fna`, `.fa`, `.fasta`) — any
other extension raises a `ValueError` before the run starts. If you supply FASTA, Tourmaline
imports both reference files itself, and `taxa_file` must be a **headerless** TSV (feature ID,
tab, taxonomy string). If you supply `.qza`, both must be artifacts.

`bt2-blca` needs plain text, so `.qza` reference files are exported to FASTA automatically for
that method.

`taxa_ranks` must match the rank structure of your reference taxonomy — it determines how the
taxonomy string is split into columns in the exported tables.

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

### REVAMP with a custom reference database (revamp-custom)

`classify_method: revamp-custom` applies REVAMP's method — BLASTn, then the lowest common
ancestor of all best hits, limited by a percent identity cutoff per rank — to your own
reference database, using **that database's taxonomy** rather than NCBI's. No NCBI
taxids, taxonkit or taxdump are involved, and the taxonomy may have any number of ranks.

```yaml
classify_method: revamp-custom
refseqs_file: /abs/path/refdb-seqs.fasta     # FASTA or .qza
taxa_file: /abs/path/refdb-taxa.tsv          # TSV (sequence ID <tab> ;-separated lineage) or .qza
taxa_ranks: domain,kingdom,phylum,subphylum,class,order,family,genus,species
revamp_dir: /abs/path/REVAMP
revamp_query_cov: 90
revamp_taxonomy_cutoffs: "97,95,90,80,70,65,60,55"   # N-1 values for N ranks, deepest first
```

**Ranks and cutoffs.** `taxa_ranks` names every field of the lineages in `taxa_file`, in
order. `revamp_taxonomy_cutoffs` takes one value per rank except the first — eight for the
nine ranks above — deepest rank first, and must not increase. As in REVAMP, an ASV whose
best hit is at or above cutoff *k* may be assigned to that rank; below the last cutoff it is
`Unassigned`, so the first rank is never assigned on its own. A wrong count stops the run
before anything executes. `revamp_blastdb` and `revamp_blast_mode` are not used.

**How it runs.** Under `{run_name}-taxonomy/revamp/`:

1. `refdb/` — `scripts/build_revamp_custom_refdb.py` gives every unique lineage an integer
   *pseudo-taxid*, renames the reference sequences `r1..rN` for `makeblastdb`, and writes
   `refdb_report.txt`: counts, rank depths, placeholder values, and **taxon names that sit
   under more than one parent** (a backbone inconsistency; ASVs hitting both sides are
   assigned above it). Read this report for any new database.
2. BLASTn against that database with REVAMP's `max_target_seqs` escalation; each hit's
   subject is replaced by its pseudo-taxid. REVAMP's `reformat_blast.R` then keeps hits
   covering `revamp_query_cov`% of the ASV and every hit at the best percent identity.
3. `scripts/revamp_lca.py` assigns taxonomy: REVAMP's assignment algorithm
   (`asv_taxonomy_processing_figureOuts.pl`) reimplemented for N ranks, since REVAMP's own
   script is fixed at seven. `ASV2Taxonomy/{run_name}_assignment_details.tsv` records each
   ASV's best percent identity, the deepest rank allowed, and its hit counts.

**Reference taxonomy handling.** Placeholder values — empty, `NA`, `Unknown`,
`unidentified`, `unassigned`, `unclassified`, or an empty prefix such as `g__` — are
treated as missing ranks, never as taxon names. Short lineages are padded with missing
trailing ranks; a lineage with more fields than `taxa_ranks` is an error. Otherwise names
are used exactly as given, rank prefixes included. As in REVAMP, a missing rank with a
named rank below it is filled (reported as `NA` in place) unless only the deepest rank is
named below it, in which case the lineage is cut there; a missing first rank makes a
reference "Environmental Unknown". When some best hits have a name at every rank down to
the allowed depth, only those are used; hits to shallower references count only when no
such hit exists.

**Differences from REVAMP's Perl.** Tested against it (`tests/test_revamp_lca.py`): on
lineages with every rank named or with internal gaps, assignments are identical. They
differ only for references missing trailing ranks, where the Perl loses track of rank
positions: it treats a reference lacking only a species as complete at species depth,
and at reduced depth removes one rank too many from such references. revamp-custom cuts
these positionally. It also does not apply REVAMP's NCBI-specific rule that `Eukaryota`
with nothing below it is "Environmental Unknown", and it compares names exactly where
REVAMP first replaces characters outside `[A-Za-z0-9 ;]` with `_`.

**BLAST run elsewhere.** Set `revamp_blast_results` to a BLASTn table produced against a
database built from `refseqs_file`, with
`-outfmt '6 qseqid pident length sseqid sacc' -subject_besthit -max_target_seqs 4000`.
Column 4 must hold `refseqs_file`'s sequence IDs (an `lcl|` prefix is accepted); the run
stops if any subject is not in the reference.

The output TSV matches `revamp`'s: third column `percent_id`, unresolved ranks `NA`, and
one rank per `taxa_ranks` entry carried through the collapsed table, barplot,
`-asv_taxa_features.tsv` and Krona plot.

Smoke test (400-sequence MiFish reference with nine ranks; BLAST takes seconds):
`./tourmaline.sh -s taxonomy -c config_03_taxonomy_revamp_custom_test.yaml -n 4`.

See [External Data](../external_data.md) for conversions if starting outside Tourmaline.

