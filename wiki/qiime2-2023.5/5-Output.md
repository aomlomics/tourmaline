## Directory structure after a completed run

A full list of outputs and associated inputs and rules for a DADA2 paired-end run in unfiltered mode, i.e., after running the command `snakemake dada2_pe_report_unfiltered`. Outputs are identical for the filtered mode except the subdirectories are named "filtered" instead of "unfiltered".

### 00-data

```
manifest_pe.csv
manifest_se.csv
metadata.tsv
repseqs_to_filter_dada2-pe.tsv
repseqs_to_filter_dada2-se.tsv
repseqs_to_filter_deblur-se.tsv
```

### 01-imported

```
fastq_counts.tsv
fastq_counts_describe.md
fastq_pe.qza
fastq_summary.qzv
refseqs.qza
reftax.qza
```

### 02-output-dada2-pe-unfiltered/00-table-repseqs

```
dada2_stats.qza
repseqs.fasta
repseqs.qza
repseqs.qzv
repseqs_amplicon_type.txt
repseqs_lengths.tsv
table.biom
table.qza
table.qzv
table_summary_features.txt
table_summary_samples.txt
```

### 02-output-dada2-pe-unfiltered/01-taxonomy

```
taxa_barplot.qzv
taxonomy.qza
taxonomy.qzv
taxonomy.tsv
```

### 02-output-dada2-pe-unfiltered/02-alignment-tree

```
aligned_repseqs.fasta
aligned_repseqs.qza
aligned_repseqs_gaps.tsv
aligned_repseqs_outliers.tsv
outliers.qza
outliers.tsv
repseqs_properties.pdf
repseqs_properties.tsv
repseqs_properties_describe.md
repseqs_to_filter_outliers.tsv
repseqs_to_filter_unassigned.tsv
rooted_tree.qza
rooted_tree.qzv
unrooted_tree.qza
```

### 02-output-dada2-pe-unfiltered/03-alpha-diversity

```
alpha_rarefaction.qzv
evenness_group_significance.qzv
evenness_vector.qza
faith_pd_group_significance.qzv
faith_pd_vector.qza
observed_features_group_significance.qzv
observed_features_vector.qza
rarefied_table.qza
shannon_group_significance.qzv
shannon_vector.qza
```

### 02-output-dada2-pe-unfiltered/04-beta-diversity

```
bray_curtis_distance_matrix.qza
bray_curtis_emperor.qzv
bray_curtis_group_significance.qzv
bray_curtis_pcoa_results.qza
jaccard_distance_matrix.qza
jaccard_emperor.qzv
jaccard_group_significance.qzv
jaccard_pcoa_results.qza
unweighted_unifrac_distance_matrix.qza
unweighted_unifrac_emperor.qzv
unweighted_unifrac_group_significance.qzv
unweighted_unifrac_pcoa_results.qza
weighted_unifrac_distance_matrix.qza
weighted_unifrac_emperor.qzv
weighted_unifrac_group_significance.qzv
weighted_unifrac_pcoa_results.qza
```

### 03-reports

```
metadata_summary.md
report_dada2-pe_unfiltered.html
report_dada2-pe_unfiltered.md
```
