#!/bin/bash

# ======= SLURM JOB CONFIGURATION =======

#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL
#SBATCH -A 
#SBATCH -p bigmem
#SBATCH -q normal
#SBATCH -N 1
#SBATCH -n 36
#SBATCH -t 48:00:00
#SBATCH -J tourmaline2_array_step3
#SBATCH --array=
#SBATCH -o /dev/null

# ======= ARGUMENT PARSING =======

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --working-directory|--wd) PROJECT_BASE="$2"; shift ;;
        --config-file) CONFIG_NAME="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# ======= VALIDATION =======

if [ -z "$PROJECT_BASE" ]; then
    echo "Error: --working-directory or --wd must be specified."
    exit 1
fi

if [ -z "$CONFIG_NAME" ]; then
    echo "Error: --config-file must be specified."
    exit 1
fi

# Convert PROJECT_BASE to absolute path
PROJECT_BASE=$(realpath "$PROJECT_BASE")

# ======= SETUP VARIABLES =======

SAMPLE_DIRS=($(ls -d "$PROJECT_BASE"/*))
THIS_SAMPLE=${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}
SAMPLE_NAME=$(basename "$THIS_SAMPLE")


# ======= CHECK FOR EXISTING OUTPUT FOLDER =======

CONFIG_FILE="${THIS_SAMPLE}/${CONFIG_NAME}"
if ! command -v yq &> /dev/null; then
    echo "yq not found. Please install yq to parse YAML files."
    exit 1
fi

# Get run name and remove any quotes
RUN_NAME=$(yq '.run_name' "$CONFIG_FILE" | tr -d '"')
EXPECTED_OUTPUT="${THIS_SAMPLE}/${RUN_NAME}_Output/${RUN_NAME}-taxonomy"

# Set the slurmout file name based on the run name
SLURMOUT_FILE="${RUN_NAME}_tourmaline2_taxonomy_${SLURM_JOB_ID}_slurmout.out"


# Redirect all output to the slurmout file
exec 1> >(tee -a "$SLURMOUT_FILE")
exec 2>&1

if [ -d "$EXPECTED_OUTPUT" ]; then
    echo "[$(date)] Skipping $SAMPLE_NAME — output already exists: $EXPECTED_OUTPUT"
else
    # ======= PIPELINE EXECUTION =======
    echo "[$(date)] Starting Taxonomy for: $SAMPLE_NAME"

    source /PATH/TO/conda.sh
    conda activate tourmaline2

    cd /PATH/TO/00_TOURMALINE_GITHUB

    # Set QIIME2 temporary directory and joblib settings
    export JOBLIB_MAX_NBYTES="5G"  # Limit memory per worker
    export OMP_NUM_THREADS=10      # Limit OpenMP threads
    export QIIME2_NUM_THREADS=10   # Limit QIIME2 threads

    # Run the pipeline with reduced number of jobs
    ./tourmaline.sh --step taxonomy --configfile "$CONFIG_FILE" --cores 10
    EXIT_STATUS=$?

    conda deactivate


    if [ $EXIT_STATUS -eq 0 ]; then
        echo "[$(date)] Finished Taxonomy for: $SAMPLE_NAME successfully."
    else
        echo "[$(date)] ERROR: Tourmaline taxonomy failed for $SAMPLE_NAME (exit code $EXIT_STATUS)."
    fi
fi


