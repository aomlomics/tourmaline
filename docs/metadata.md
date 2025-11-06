## Bioinformatics Metadata

Generate a single TSV with metadata about your analysis using FAIR eDNA terms.

### Usage

```bash
python scripts/format_analysisMetadata.py \
  -s config_01_qaqc.yaml \
  -r config_02_repseqs.yaml \
  -t config_03_taxonomy.yaml \
  -p my_project \
  -o my-tourmaline-metadata.tsv
``;

Options (subset):

```text
-a, --assay_name               override assay_name (else from samples config)
-A, --analysis_run_name        override analysis_run_name (else from samples config)
-T, --tourmaline_metadata      path to tourmaline metadata (if running outside repo)
```

See `--help` for the full CLI. Ensure the three config files match the runs whose outputs you want to document.


