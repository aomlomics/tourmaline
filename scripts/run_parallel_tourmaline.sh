#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 --config-dir DIR --config-prefix PREFIX --step STEP --parallel-jobs N --cores-per-job M"
    echo "Example: $0 --config-dir parameter_sweep_configs --config-prefix config-01-qaqc --step qaqc --parallel-jobs 4 --cores-per-job 6"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --config-dir) config_dir="$2"; shift ;;
        --config-prefix) config_prefix="$2"; shift ;;
        --step) step="$2"; shift ;;
        --parallel-jobs) parallel_jobs="$2"; shift ;;
        --cores-per-job) cores_per_job="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if required parameters are provided
if [ -z "$config_dir" ] || [ -z "$config_prefix" ] || [ -z "$step" ] || [ -z "$parallel_jobs" ] || [ -z "$cores_per_job" ]; then
    usage
fi

# Check if GNU parallel is installed
if ! command -v parallel &> /dev/null; then
    echo "GNU parallel is required but not installed. Please install it first."
    exit 1
fi

# Check if the config directory exists
if [ ! -d "$config_dir" ]; then
    echo "Config directory $config_dir does not exist"
    exit 1
fi

# Function to run tourmaline for a single config
run_tourmaline() {
    config_file="$1"
    step="$2"
    cores="$3"
    
    echo "Processing $config_file..."
    ./tourmaline.sh --step "$step" --configfile "$config_file" --cores "$cores"
}
export -f run_tourmaline

# Find all config files matching the prefix and run them in parallel
find "$config_dir" -name "${config_prefix}_*.yaml" | \
    parallel --jobs "$parallel_jobs" run_tourmaline {} "$step" "$cores_per_job"

echo "All jobs completed!"
