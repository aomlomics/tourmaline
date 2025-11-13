## Taxonomy Step

Assigns taxonomy and generates visualizations using one of four methods:

- Naive Bayes (QIIME 2 classify-sklearn)
- Consensus BLAST
- Consensus VSEARCH
- BT2-BLCA (Anacapa)

### Configuration

```yaml
classify_method: naive-bayes        # or consensus-blast, consensus-vsearch, bt2-blca
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

# Optional metadata for barplots
sample_metadata_file: 00-data/metadata.tsv
```

### Outputs

- Taxonomy assignments and visualizations (barplots, summaries)
- Optional collapsed table at `collapse_taxalevel`
- Combined ASV/taxonomy/sequence TSV for downstream FAIR metadata workflows
- Exported taxonomy TSV alongside the QIIME 2 taxonomy artifact

See [External Data](../external_data.md) for conversions if starting outside Tourmaline.

