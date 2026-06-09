Example usage of the updated format_analysisMetadata.py script.

The script now takes a working directory and run names for each step,
then automatically finds the config files within each step's folder.
It also copies additional output files to an output folder with the analysis_run_name prefix.

Example directory structure:
working_dir/
├── my-run-qaqc/
│   └── my-run-qaqc_config.yaml
├── my-run-repseqs/
│   ├── my-run-repseqs_config.yaml
│   └── my-run-table.tsv
└── my-run-taxonomy/
    ├── my-run-taxonomy_config.yaml
    └── my-run-asv_taxa_features.tsv

Output folder will contain:
├── my-analysis_metadata.tsv
├── my-analysis_asv_taxa_features.tsv
└── my-analysis_table.tsv

Usage:
python format_analysisMetadata.py \
    --working_dir /path/to/working_dir \
    --qaqc_run_name my-run \
    --repseqs_run_name my-run \
    --taxonomy_run_name my-run \
    --project_id my-project \
    --output_folder /path/to/output