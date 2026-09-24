## Quick Start

Tourmaline 2 provides a modular workflow for processing amplicon sequencing data. The pipeline consists of three main steps that can be run together or independently:

1. **QA/QC** - Process raw FASTQ files, optional primer trimming, generate QIIME 2 artifact
2. **Repseqs** - Generate ASVs using DADA2 or Deblur, optional filtering, produce feature table and representative sequences
3. **Taxonomy** - Assign taxonomy using one of five methods, generate visualizations

### Getting Started

1. **Install and Setup**: See [Install and Setup](install.md) for requirements and environment setup.
   - Install QIIME 2 (2024.10) amplicon workflow
   - Create Snakemake conda environment
   - Clone Tourmaline repository

2. **Configuration**: See [Configuration](configuration.md) for setting up config files for each step.
   - Create config files for each step you plan to run
   - Configure parameters for your data type and analysis needs

3. **Running**: See [Running](running.md) for how to use `tourmaline.sh` and examples.
   - Run single steps or all steps together
   - Use the `tourmaline.sh` script

### Example Workflow

```bash
# Activate environment
conda activate snakemake-tour2

# Run all steps
./tourmaline.sh -s qaqc,repseqs,taxonomy \
  -c config_01_qaqc.yaml,config_02_repseqs.yaml,config_03_taxonomy.yaml \
  -n 6

# Or run a single step
./tourmaline.sh -s taxonomy -c config_03_taxonomy.yaml -n 6
```

## Directory structure

The pipeline creates the following directory structure for outputs:

```
output_dir/
├── [run_name]-qaqc/       # QA/QC outputs
├── [run_name]-repseqs/    # Representative sequences outputs
└── [run_name]-taxonomy/   # Taxonomy assignment outputs
```

Each directory contains the relevant outputs for that step of the pipeline, along with a copy
of the config file that produced them (`{run_name}-{step}_config.yaml`).

### Before you start

Tourmaline expects **demultiplexed** FASTQ files — one file per sample, or two for paired-end
data. You will also need your primer sequences (if you want them trimmed) and a reference
database for your marker gene.

Tips for a first run:

- Try the small example data in `00-data/` first, to confirm your install works.
- Run one step at a time and inspect the `.qzv` outputs at [view.qiime2.org](https://view.qiime2.org)
  before moving on. Look at the Step 1 quality plots *before* choosing DADA2 truncation lengths.
- Add `--dryrun` to a direct Snakemake call to see what would run without running it.
- Give each attempt its own `run_name`, so runs sit side by side and can be compared.

### Next steps

- [Analysis Metadata](metadata.md): produce FAIR eDNA metadata for the NOAA Ocean DNA Explorer
- [Troubleshooting](troubleshooting.md): when something goes wrong

