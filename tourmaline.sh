#!/bin/bash
##
# @Description: Runs the snakemake pipeline Tourmaline 2.0, either all at once or modularly.
# Requires snakemake env activated. 
#command: ./tourmaline.sh -s taxonomy --configfile config-03-taxonomy.yaml --cores 6
# test full: ./tourmaline.sh -s qaqc,repseqs,taxonomy --configfile config-01-sample.yaml,config-02-repseqs.yaml,config-03-taxonomy.yaml --cores 6
##



# Function to display usage
usage() {
    echo "Usage: $0 --step qaqc,repseqs,taxonomy --configfile config1,config2,... --cores N"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --step|-s) steps="$2"; shift ;;
        --configfile|-c) configfiles="$2"; shift ;;
        --cores|-n) cores="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if required parameters are provided
if [ -z "$steps" ] || [ -z "$configfiles" ] || [ -z "$cores" ]; then
    usage
fi

# Split the steps and configfiles into arrays
IFS=',' read -r -a step_array <<< "$steps"
IFS=',' read -r -a configfile_array <<< "$configfiles"

# Check if the number of steps matches the number of config files
if [ "${#step_array[@]}" -ne "${#configfile_array[@]}" ]; then
    echo "The number of steps must match the number of config files."
    usage
fi

# Iterate over each step and run the corresponding code chunk
for index in "${!step_array[@]}"; do
    step="${step_array[$index]}"
    CONFIG="${configfile_array[$index]}"
    case $step in
        qaqc)
            echo "Running QA/QC step with configfile $configfile and cores $cores\n"
            paired=$(yq -r '.paired_end' $CONFIG);
            trim=$(yq -r '.to_trim' $CONFIG);
            if [[ "${paired}" = true ]]; then
                if [[ "${trim}" = true ]]; then
                #snakemake --use-conda -s qaqc_step_paired.Snakefile --configfile $CONFIG --cores $CORES  --latency-wait 15
                    snakemake --use-conda -s qaqc_step_paired.Snakefile trim_pe_all --configfile $CONFIG --cores $cores  --latency-wait 15
                else
                    snakemake --use-conda -s qaqc_step_paired.Snakefile no_trim_pe_all --configfile $CONFIG --cores $cores --latency-wait 15
                fi;
            else
                #snakemake --use-conda -s qaqc_step_single.Snakefile --configfile $CONFIG --cores $CORES --latency-wait 15
                if [[ "${trim}" = true ]]; then
                    snakemake --use-conda -s qaqc_step_single.Snakefile trim_se_all --configfile $CONFIG --cores $cores  --latency-wait 15
                else
                    snakemake --use-conda -s qaqc_step_single.Snakefile no_trim_se_all --configfile $CONFIG --cores $cores  --latency-wait 15
                fi;
            fi;
            ;;
        repseqs)
            echo "Running repseqs step with configfile $CONFIG and cores $cores\n"
            asv_method=$(yq -r '.asv_method' $CONFIG);
            if [[ ${asv_method} = dada2pe ]]; then
                snakemake --use-conda -s repseqs_step.Snakefile run_dada2_pe_denoise --configfile $CONFIG --cores $cores  --latency-wait 15
            elif [[ "${asv_method}" = "dada2se" ]]; then
                snakemake --use-conda -s repseqs_step.Snakefile run_dada2_se_denoise --configfile $CONFIG --cores $cores  --latency-wait 15
            elif [[ "${asv_method}" = "deblur" ]]; then
                snakemake --use-conda -s repseqs_step.Snakefile run_deblur_se_denoise --configfile $CONFIG --cores $cores  --latency-wait 15
            else
                echo "ASV method not recognized"
            fi;
            ;;
        taxonomy)
            echo "Running taxonomy step with configfile $CONFIG and cores $cores"
            echo -n "Running taxonomy step\n"
            snakemake --use-conda -s taxonomy_step.Snakefile --configfile $CONFIG --cores $cores  --latency-wait 15
            ;;
        *)
            echo "Unknown step: $step"
            usage
            ;;
    esac
done



