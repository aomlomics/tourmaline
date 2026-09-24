## Bioinformatics Metadata

`scripts/format_analysisMetadata.py` generates a single TSV describing your analysis using
[FAIR eDNA](https://fair-edna.github.io/) terms. The result can be read into the
[NOAA Ocean DNA Explorer](https://www.ngi.msstate.edu/node).

This is a standalone script, not a `tourmaline.sh` step. Run it after the qaqc, repseqs and
taxonomy steps have finished.

### How it finds your run

You give it the **results directory** and the **run name** of each step — not config file
paths. Each Tourmaline step copies its config into its own output directory as
`{run_name}-{step}_config.yaml`, and the script reads those copies, so the metadata always
describes the parameters a run actually used.

It expects the standard output layout:

```
working_dir/
├── my-qaqc-run-qaqc/
│   └── my-qaqc-run-qaqc_config.yaml
├── my-repseqs-run-repseqs/
│   ├── my-repseqs-run-repseqs_config.yaml
│   └── my-repseqs-run-table.tsv
└── my-tax-run-taxonomy/
    ├── my-tax-run-taxonomy_config.yaml
    └── my-tax-run-asv_taxa_features.tsv
```

The three run names may differ — that is the normal case when you reuse one QA/QC run across
several denoising or classification runs.

### Usage

Run from the Tourmaline directory, so the default `-T ./00-data/tourmaline_metadata.yaml`
resolves:

```bash
python scripts/format_analysisMetadata.py \
  -w ../v2-results \
  -q my-qaqc-run \
  -r my-repseqs-run \
  -t my-tax-run \
  -p my_project \
  -O output_folder/
```

### Options

```text
-w, --working_dir          (required) directory holding the step output folders
-q, --qaqc_run_name        (required) run name of the qaqc step
-r, --repseqs_run_name     (required) run name of the repseqs step
-t, --taxonomy_run_name    (required) run name of the taxonomy step
-p, --project_id           (required) value for project_id
-O, --output_folder        (required) folder where outputs are written
-a, --assay_name           override assay_name (else taken from the qaqc config)
-A, --analysis_run_name    override analysis_run_name (else the taxonomy run name)
-T, --tourmaline_metadata  path to the tourmaline metadata YAML
                           (default: ./00-data/tourmaline_metadata.yaml)
```

Run with `--help` for the authoritative list.

### Outputs

The output folder receives three files, each prefixed with the analysis run name:

```
output_folder/
├── {analysis_run_name}_analysisMetadata.tsv   # the FAIR eDNA metadata table
├── {analysis_run_name}_asv_taxa_features.tsv  # copied from the taxonomy run
└── {analysis_run_name}_table.tsv              # copied from the repseqs run
```

### Notes

- `assay_name` comes from the qaqc config unless overridden. Use a term from the
  [NOAA Omics metabarcoding assays](https://github.com/NOAA-Omics/noaa-omics-metabarcoding-assays/blob/main/assays.tsv)
  controlled vocabulary; open an
  [issue](https://github.com/NOAA-Omics/noaa-omics-metabarcoding-assays/issues) if yours is missing.
- Software versions are read from `00-data/tourmaline_metadata.yaml`. If you upgrade QIIME 2 or
  a classifier, update that file so the recorded versions stay accurate.
- The script stops if a step directory or its config cannot be found — check that the run
  names match the directory names under `--working_dir`. A missing `asv_taxa_features.tsv` or
  `-table.tsv` only prints a warning: the metadata TSV is still written, without that copy.
