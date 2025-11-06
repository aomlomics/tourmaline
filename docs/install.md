## Install and Setup

Tourmaline 2 uses [QIIME 2 2024.10](https://docs.qiime2.org/2024.10/install/) for analysis and [Snakemake](https://snakemake.readthedocs.io/en/stable/) for workflow management.

### Requirements

- Conda (Miniconda recommended)
- QIIME 2 (2024.10) amplicon workflow: environment name `qiime2-amplicon-2024.10`
- Snakemake conda environment (name example: `snakemake-tour2`)

### Create environments

```bash
conda create -c conda-forge -c bioconda -n snakemake-tour2 snakemake biopython yq parallel
```

Install QIIME 2 (see QIIME 2 docs for OS specifics), ensuring the environment is named exactly:

```bash
conda create -n qiime2-amplicon-2024.10  # or use official installer for 2024.10
```

Optional (BLCA method only):

```bash
conda create -c conda-forge -c bioconda -n bt2-blca biopython muscle=3.8 bowtie2
```

### Get Tourmaline 2

```bash
git clone https://github.com/aomlomics/tourmaline.git
cd tourmaline
```

### Running requirements

- Activate `snakemake-tour2` when running the pipeline
- Keep `qiime2-amplicon-2024.10` installed (used by Snakemake rules via `--use-conda`)

### Legacy v1

For the legacy pipeline and docs, see the V1 branch and legacy docs referenced in [Citation & Legacy](citation_legacy.md).


