## Install and Setup

Tourmaline 2 uses [QIIME 2 2024.10](https://docs.qiime2.org/2024.10/install/) for analysis and
[Snakemake](https://snakemake.readthedocs.io/en/stable/) for workflow management.

Snakemake rules request conda environments **by name**, so the names below must match exactly —
they are not created for you from environment files.

### Requirements

- Conda ([Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) recommended)
- QIIME 2 2024.10 amplicon distribution, in an environment named `qiime2-amplicon-2024.10`
- A Snakemake environment, named `snakemake-tour2`

### Create environments

**QIIME 2** — follow the [official 2024.10 install instructions](https://docs.qiime2.org/2024.10/install/)
for your operating system. The default environment name from the official installer is already
`qiime2-amplicon-2024.10`; if you name it something else, the rules will not find it. You do not
need to install anything else into this environment (except for the tax-credit step, below).

Verify:

```bash
conda activate qiime2-amplicon-2024.10
qiime --version
```

**Snakemake**:

```bash
conda create -c conda-forge -c bioconda -n snakemake-tour2 snakemake biopython yq parallel
```

This is the environment you activate to *run* Tourmaline. It includes GNU `parallel`, used by the
parameter sweep tooling.

### Optional environments, per feature

Create these only if you use the corresponding feature.

**BLCA taxonomy** (`classify_method: bt2-blca`):

```bash
conda create -c conda-forge -c bioconda -n bt2-blca biopython muscle=3.8 bowtie2
```

**REVAMP taxonomy** (`classify_method: revamp`):

```bash
conda create -c conda-forge -c bioconda -n revamp \
  "blast>=2.13" "taxonkit>=0.20" r-base r-dplyr bioconductor-biostrings \
  perl perl-list-moreutils krona
```

`taxonkit` must be 0.20 or newer. REVAMP also needs a clone of
[REVAMP](https://github.com/McAllister-NOAA/REVAMP) and a local NCBI `nt` database with prepared
taxonomy files — see [Taxonomy step](steps/taxonomy.md#revamp).

**Krona plots** (`make_krona: True`, any classify method):

```bash
conda create -c conda-forge -c bioconda -n krona krona
```

No Krona taxonomy database download is required.

**Tax-credit benchmarking** (the `tax-credit` step) needs the sibling
[tax-credit](https://github.com/aomlomics/tax-credit) package installed into the QIIME 2
environment:

```bash
conda activate qiime2-amplicon-2024.10
pip install -e ../tax-credit
```

### Get Tourmaline 2

The default branch is `V2`:

```bash
git clone https://github.com/aomlomics/tourmaline.git
cd tourmaline
```

The `tax-credit` step currently lives on the `feature/tax-credit-module` branch, which is not yet
merged into `V2`.

### Running requirements

- Activate `snakemake-tour2` before running anything.
- Run from the Tourmaline repository root — rules invoke `scripts/...` by relative path. Outputs
  can go anywhere, via `output_dir` in the config.
- Keep `qiime2-amplicon-2024.10` installed; rules pull it via `--use-conda`.

### Check your install

Example data in `00-data/` lets you confirm the install before pointing at real data:

```bash
conda activate snakemake-tour2
./tourmaline.sh -s qaqc -c config_01_qaqc.yaml -n 4
```

Add `--dryrun` to a direct Snakemake call to validate a config without running anything — see
[Running](running.md).

### Legacy v1

For the legacy pipeline and its docs, see the
[V1 branch](https://github.com/aomlomics/tourmaline/tree/V1) and
[Citation & Legacy](citation_legacy.md).
