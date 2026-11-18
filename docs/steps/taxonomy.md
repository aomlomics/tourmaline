## Taxonomy Step

Assigns taxonomy and generates visualizations using one of six methods:

- Naive Bayes (QIIME 2 classify-sklearn)
- Consensus BLAST
- Consensus VSEARCH
- BT2-BLCA (Anacapa)
- REVAMP (BLASTn against a local NCBI nt database; see [REVAMP](#revamp) below)
- REVAMP with your own reference database and taxonomy, any number of ranks
  (`revamp-custom`; see [below](#revamp-with-your-own-reference-database-revamp-custom))

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

# REVAMP options (classify_method: revamp or revamp-custom)
revamp_dir: /abs/path/REVAMP              # REVAMP clone (both)
revamp_blastdb: /abs/path/blastdb         # revamp only: nt volumes plus a prepared taxdump/
revamp_blast_results: /abs/path/ASV_blastn_nt.btab  # optional: BLAST run elsewhere
revamp_blast_mode: mostEnvOUT             # revamp only: allIN | allEnvOUT | mostEnvOUT
revamp_query_cov: 90                      # percent of ASV length a hit must cover
revamp_taxonomy_cutoffs: "97,95,90,80,70,60"  # percent ID cutoffs, deepest rank first:
                                          # revamp: S,G,F,O,C,P; revamp-custom: one per rank but the first

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

### REVAMP with your own reference database (revamp-custom)

`classify_method: revamp-custom` assigns taxonomy with REVAMP's method — BLASTn, then the
lowest common ancestor of all best hits, limited by a percent identity cutoff for each
rank — against **your own reference database, using its own taxonomy**. Use it when you
have a curated reference (a regional 12S database, PR2, MIDORI, ...) and want results
named in that database's taxonomy rather than NCBI's.

Compared with `revamp`:

| | `revamp` | `revamp-custom` |
|---|---|---|
| Reference sequences | local NCBI `nt` | `refseqs_file` |
| Taxonomy | NCBI, via taxids and taxonkit | `taxa_file`, as written |
| Ranks | always the seven kingdom..species | any number, from `taxa_ranks` |
| Cutoffs (`revamp_taxonomy_cutoffs`) | six | one fewer than the number of ranks |
| Needs `revamp_blastdb` / taxdump | yes | no |
| Environmental filtering (`revamp_blast_mode`) | yes | not applicable |

#### 1. What you need

- The `revamp` conda environment and a clone of
  [REVAMP](https://github.com/McAllister-NOAA/REVAMP) — see [Install](../install.md). No
  NCBI `nt` database or taxonomy dump is needed.
- Your reference sequences (`refseqs_file`): FASTA or a QIIME 2 `FeatureData[Sequence]`
  `.qza`.
- Their taxonomy (`taxa_file`): a tab-separated file or a `FeatureData[Taxonomy]` `.qza`.

#### 2. Prepare the reference files

The taxonomy file has one line per reference sequence: the sequence ID, a tab, then the
lineage with ranks separated by `;`, highest rank first. A header line (`Feature ID	Taxon`)
and extra columns (such as `Confidence`) are fine.

```
Feature ID	Taxon
AB282828.1	Eukaryota;Metazoa;Chordata;Vertebrata;Actinopteri;Lophiiformes;Antennariidae;Antennarius;Antennarius striatus
AJM0165	Eukaryota;Metazoa;Chordata;Vertebrata;Actinopteri;NA;Lutjanidae;Lutjanus;Lutjanus analis
```

- **Every lineage should have one field per rank**, in the same order. A lineage with
  fewer fields is padded with missing ranks at the end; one with more fields than
  `taxa_ranks` stops the run.
- **Mark a missing rank with an empty field or a placeholder.** Empty fields, `NA`, `N/A`,
  `Unknown`, `unidentified`, `unassigned`, `unclassified` and an empty rank prefix such as
  `g__` all count as missing — never as a taxon name, so two unrelated `NA` genera are not
  treated as the same genus. Keep the field so later ranks stay in position:
  `...;Actinopteri;NA;Lutjanidae;...`, not `...;Actinopteri;Lutjanidae;...`.
- **Names are used exactly as written**, rank prefixes (`g__Antennarius`) included, and
  come back in the output the same way.
- Sequence IDs must match between the two files; sequences with no taxonomy are left out
  of the database and counted in the report.

#### 3. Configure

```yaml
classify_method: revamp-custom
database_name: mifish-gomex-12S                 # used in the analysis metadata
refseqs_file: /abs/path/refdb-seqs.fasta
taxa_file: /abs/path/refdb-taxa.tsv
taxa_ranks: domain,kingdom,phylum,subphylum,class,order,family,genus,species
collapse_taxalevel: 9                           # 1 .. number of taxa_ranks
revamp_dir: /abs/path/REVAMP
revamp_query_cov: 90
revamp_taxonomy_cutoffs: "97,95,90,80,70,65,60,55"
```

`taxa_ranks` names every rank in your taxonomy, in order. `revamp_taxonomy_cutoffs` gives
one percent identity per rank **except the first**, starting from the deepest rank — so
nine ranks take eight values. For the example above:

| Best-hit percent identity | Deepest rank the ASV can be assigned |
|---|---|
| ≥ 97 | species |
| 95 – < 97 | genus |
| 90 – < 95 | family |
| 80 – < 90 | order |
| 70 – < 80 | class |
| 65 – < 70 | subphylum |
| 60 – < 65 | phylum |
| 55 – < 60 | kingdom |
| < 55 | `Unassigned` (the first rank is never assigned on its own, as in REVAMP) |

The cutoffs only set the *deepest possible* rank: when the best hits disagree, the
assignment stops at the rank they share. Values must not increase from left to right. A
wrong count stops the run before anything executes, with a message saying how many are
needed. `revamp_query_cov` (0-100) is the share of the ASV's length a hit must cover to
count. `revamp_blastdb` and `revamp_blast_mode` are ignored.

#### 4. Run

```bash
snakemake --use-conda -s taxonomy_step.Snakefile --configfile config_03_taxonomy.yaml --cores 4 --dryrun
./tourmaline.sh -s taxonomy -c config_03_taxonomy.yaml -n 4
```

BLAST uses `classify_threads`. The BLAST database is built from `refseqs_file` on each
run; for the 400-sequence test reference, building and searching take seconds.

#### 5. Check the reference report

Read `{run_name}-taxonomy/revamp/refdb/refdb_report.txt` the first time you use a
database. It lists the sequence and lineage counts, the deepest named rank per sequence,
which placeholder values were treated as missing, and warnings such as:

```
WARNINGS:
  - 6 taxon names occur under more than one parent lineage. ASVs hitting more than one of
    them are assigned above the conflicting rank. Examples:
      family 'Pomacanthidae' under 2 parents, e.g. ...;Actinopteri; | ...;Actinopteri;Acanthuriformes
```

A name under two parents is an inconsistency in the reference taxonomy — here one record
lacks the order. Wherever the best hits include both versions, the assignment stops above
the conflict, so fixing these in `taxa_file` gives deeper, more consistent assignments.
The report also flags sequences with no first rank (treated as "Environmental Unknown",
as in REVAMP) and gaps that cannot be filled (for example a missing genus with a species
present, which cuts that lineage at the gap).

#### 6. Read the results

The outputs are the same as for every method: `{run_name}-taxonomy.tsv` / `.qza`, the
barplot, the collapsed table, `-asv_taxa_features.tsv` and the optional Krona plot — with
**one rank per `taxa_ranks` entry**. As with `revamp`:

- The third column is `percent_id`, the best hit's percent identity (0-100), not a 0-1
  `Confidence`.
- A rank that could not be resolved, or that is missing in the reference, is `NA` in place;
  trailing `NA`s are dropped, so a genus never shifts into the family slot.
- ASVs with no hit covering `revamp_query_cov`% of their length, or below the last cutoff,
  are `Unassigned`.

`revamp/ASV2Taxonomy/{run_name}_assignment_details.tsv` explains each assignment: the best
percent identity, the deepest rank the cutoffs allowed, the number of best hits, and how
many of them had a full lineage down to that rank. Compare `max_depth_rank` with the
assigned taxonomy to see whether an ASV stopped early because of its identity or because
its best hits disagreed.

#### Choosing cutoffs

REVAMP suggests `97,95,90,80,70,60` for rRNA genes and `95,92,87,77,67,60` for
protein-coding genes (species to phylum); for a taxonomy with more ranks, place values for
the extra ranks between those of their neighbours. The best values depend on the marker and
the database, so benchmark them with the [tax-credit step](tax_credit.md), which can sweep
several cutoff sets and supports `revamp-custom` for every evaluation method.

#### BLAST run elsewhere

To BLAST on another machine, build a database from `refseqs_file` and run

```bash
blastn -db <refdb> -query ASVs.fa -outfmt '6 qseqid pident length sseqid sacc' \
  -subject_besthit -max_target_seqs 4000 -out ASV_blastn_custom.btab
```

then set `revamp_blast_results` to the output. Query IDs must be the feature IDs of your
representative sequences, and column 4 must hold `refseqs_file`'s sequence IDs (an `lcl|`
prefix is accepted). The run stops if any subject is not in the reference.

#### Things to know

- **Distant ASVs may get no hit.** BLASTn runs in its default megablast mode, as in REVAMP,
  which is tuned for close matches. ASVs from groups poorly represented in the reference
  can find no hit covering `revamp_query_cov`% of their length and come back `Unassigned`
  rather than assigned at a high rank.
- **Results are only as good as the reference taxonomy.** Inconsistent lineages (see the
  report) and missing ranks limit how deep assignments can go.

#### How it works

Under `{run_name}-taxonomy/revamp/`:

1. `refdb/` — `scripts/build_revamp_custom_refdb.py` gives every unique lineage an integer
   *pseudo-taxid* (standing in for an NCBI taxid), renames the sequences `r1..rN` for
   `makeblastdb` (written as local IDs, `lcl|r1`, so BLAST never mistakes one for a
   GenBank accession), and writes the report.
2. BLASTn against that database, with REVAMP's escalation of `-max_target_seqs` until
   every best hit fits; each hit's subject becomes its pseudo-taxid. REVAMP's
   `reformat_blast.R` then keeps hits covering `revamp_query_cov`% of the ASV and every
   hit at the best percent identity.
3. `scripts/revamp_lca.py` assigns taxonomy with REVAMP's algorithm (from
   `asv_taxonomy_processing_figureOuts.pl`), rewritten for any number of ranks because
   REVAMP's script is fixed at seven. As in REVAMP, a missing rank with a named rank below
   it is filled, and when some best hits have a name at every rank down to the allowed
   depth, only those are combined.

`tests/test_revamp_lca.py` checks it against REVAMP's own script on seven-rank input: on
lineages with every rank named, or with internal gaps, the assignments are identical. They
differ only for references missing their last ranks, where REVAMP's script loses track of
rank positions (it treats a reference lacking only a species as complete at species
depth, and at reduced depth drops one rank too many); revamp-custom cuts these by
position. It also skips REVAMP's NCBI-specific rule that `Eukaryota` with nothing below it
is "Environmental Unknown", and compares names exactly where REVAMP first replaces
characters outside `[A-Za-z0-9 ;]` with `_`.

Smoke test (400-sequence MiFish reference with nine ranks; BLAST takes seconds):
`./tourmaline.sh -s taxonomy -c config_03_taxonomy_revamp_custom_test.yaml -n 4`.

See [External Data](../external_data.md) for conversions if starting outside Tourmaline.

