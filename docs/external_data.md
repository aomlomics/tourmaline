## External Data

Tourmaline 2 can start steps from externally generated inputs.

### Repseqs step: external demultiplexed sequences

Provide a QIIME 2 artifact of demultiplexed sequences:

```yaml
# in config_02_repseqs.yaml
fastq_qza_file: /abs/path/to/fastq.qza
```

Or reuse QA/QC outputs by specifying `sample_run_name`.

### Taxonomy step: external repseqs and table

Provide both files:

```yaml
# in config_03_taxonomy.yaml
repseqs_qza_file: /abs/path/to/repseqs.qza
table_qza_file: /abs/path/to/table.qza
```

### Conversions

Convert FASTA to QZA (ASV sequences):

```bash
qiime tools import \
   --type 'FeatureData[Sequence]' \
   --input-path my-asvs.fasta \
   --output-path output-asvs.qza
```

Convert BIOM v2.1 to QZA (feature table):

```bash
qiime tools import \
  --input-path feature-table-v210.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path feature-table.qza
```

If starting from a TSV count table, first convert to BIOM (see BIOM docs), then import.


