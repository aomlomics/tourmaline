## External Data

Unlike Tourmaline 1, **any step can be the starting point**. Steps chain through files on disk, so
a step will happily read an artifact produced by another program — as long as it is a QIIME 2
artifact of the type that step expects.

| You already have | Start at | Config keys |
|---|---|---|
| Demultiplexed FASTQs, untrimmed or trimmed | qaqc | `raw_fastq_path`, `trimmed_fastq_path`, `sample_manifest_file` |
| An imported demultiplexed `.qza` | qaqc or repseqs | `preexisting_fastq_qza` (qaqc) or `fastq_qza_file` (repseqs) |
| ASV sequences and a count table | taxonomy | `repseqs_qza_file` **and** `table_qza_file` |

Run the conversions below in the QIIME 2 environment:

```bash
conda activate qiime2-amplicon-2024.10
```

> Artifacts written by one QIIME 2 version are often **not** readable by earlier versions. Import
> with 2024.10 if the result will be used by Tourmaline 2.

### QA/QC step: already-imported sequences

If your reads are already a QIIME 2 demultiplexed artifact, skip importing:

```yaml
# in config_01_qaqc.yaml
preexisting_fastq_qza: /abs/path/to/demux.qza
```

### Repseqs step: external demultiplexed sequences

```yaml
# in config_02_repseqs.yaml
fastq_qza_file: /abs/path/to/fastq.qza
```

Or reuse QA/QC outputs by specifying `qaqc_run_name`.

### Taxonomy step: external repseqs and table

Provide **both**:

```yaml
# in config_03_taxonomy.yaml
repseqs_qza_file: /abs/path/to/repseqs.qza
table_qza_file: /abs/path/to/table.qza
```

The two must describe the same set of ASV IDs.

### Conversions

#### FASTQ files → demultiplexed `.qza`

Needs a manifest file mapping sample names to absolute file paths — see
[QA/QC step](steps/qaqc.md#manifest-formats) for the format.

Paired-end:

```bash
qiime tools import \
   --type 'SampleData[PairedEndSequencesWithQuality]' \
   --input-path my_pe.manifest \
   --output-path output-file_pe_fastq.qza \
   --input-format PairedEndFastqManifestPhred33V2
```

Single-end:

```bash
qiime tools import \
   --type 'SampleData[SequencesWithQuality]' \
   --input-path my_se.manifest \
   --output-path output-file_se_fastq.qza \
   --input-format SingleEndFastqManifestPhred33V2
```

#### FASTA → `.qza` (ASV sequences)

```bash
qiime tools import \
   --type 'FeatureData[Sequence]' \
   --input-path my-asvs.fasta \
   --output-path output-asvs.qza
```

#### BIOM → `.qza` (feature table)

[Check the BIOM version first](https://docs.qiime2.org/2024.10/tutorials/importing/#feature-table-data),
then use the matching `--input-format`:

```bash
# BIOM v2.1 (HDF5)
qiime tools import \
  --input-path feature-table-v210.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path feature-table.qza

# BIOM v1.0.0 (JSON)
qiime tools import \
  --input-path feature-table-v100.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path feature-table.qza
```

#### TSV count table → `.qza`

For a TSV with unique sequences (or feature IDs) as rows and samples as columns, convert to BIOM
first:

```bash
biom convert -i otu_table.txt -o new_otu_table.biom --to-hdf5 --table-type="OTU table"

qiime tools import \
  --input-path new_otu_table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path feature-table.qza
```

#### Reference databases

`refseqs_file` accepts either a QIIME 2 artifact (`.qza`) or plain FASTA (`.fna`, `.fa`,
`.fasta`); any other extension stops the run with a `ValueError`. When you supply FASTA,
Tourmaline imports both files for you, and `taxa_file` must then be a **headerless** TSV
(`HeaderlessTSVTaxonomyFormat`) — feature ID and taxonomy string, tab separated, no header row.

For `bt2-blca` the conversion runs the other way: that method needs FASTA and plain text, so
`.qza` inputs are exported automatically.

To import them yourself instead:

```bash
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ref_seqs.fasta \
  --output-path ref_seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path ref_taxa.tsv \
  --output-path ref_taxa.qza
```

Make sure `taxa_ranks` in your taxonomy config matches the rank structure of the taxonomy file.

### Checking an artifact

```bash
qiime tools peek my-artifact.qza     # type and UUID
qiime tools export --input-path my-artifact.qza --output-path exported/
```
