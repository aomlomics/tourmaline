#!/bin/bash
#SBATCH --job-name=tour2-tax-credit
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=tourmaline_tax_credit_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source "${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"
conda activate snakemake-tour2
./tourmaline.sh --step tax-credit --configfile config_04_tax_credit.yaml --cores "${SLURM_CPUS_PER_TASK:-8}"
