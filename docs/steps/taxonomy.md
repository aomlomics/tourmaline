## Taxonomy Step

Assigns taxonomy and generates visualizations using one of four methods:

- Naive Bayes (QIIME 2 classify-sklearn)
- Consensus BLAST
- Consensus VSEARCH
- BT2-BLCA (Anacapa)

### Configuration

```yaml
classify_method: naive-bayes  # or consensus-blast, consensus-vsearch, bt2-blca
classify_threads: 4
collapse_taxalevel: 0

# Inputs from Repseqs step or external
repseqs_run_name: repseqs_run
repseqs_qza_file: /abs/path/to/repseqs.qza
table_qza_file: /abs/path/to/table.qza

# Reference database
database_name: PR2
refseqs_file: 00-data/refseqs.fna
taxa_file: 00-data/reftax.tsv
pretrained_classifier: /abs/path/to/classifier.qza  # optional for naive-bayes
bowtie_database: /abs/path/to/bt2/index/            # only for bt2-blca
taxa_ranks: kingdom,phylum,class,order,family,genus,species

# Thresholds (examples)
skl_confidence: 0.7
perc_identity: 0.8
query_cov: 0.8
min_consensus: 0.51
confidence_thres: 0.8  # bt2-blca
```

### Outputs

- Taxonomy assignments and visualizations (barplots, summaries)
- Optional collapsed table at `collapse_taxalevel`

See [External Data](../external_data.md) for conversions if starting outside Tourmaline.

