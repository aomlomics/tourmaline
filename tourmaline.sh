#!/bin/bash
##
# @Description: Runs the snakemake pipeline Tourmaline 2.0, either all at once or modularly.
#command: ./tourmaline.sh -s sample --configfile config-01-sample.yaml --cores 4
##

function show_usage (){
    printf "Usage: $0 [options [parameters]]\n"
    printf "\n"
    printf "Options:\n"
    printf " -s|--step, sample\n"
    printf " -c|--configfile, config-01-sample.yaml\n"
    printf " -n|--cores, numeric\n"
    printf " -h|--help, Print help\n"

return 0
}


while [ ! -z "$1" ]; do
  case "$1" in
     --step|-s)
         shift
         STEP="$1"
         ;;
     --configfile|-c)
         shift
         CONFIG="$1"
         ;;
     --cores|-n)
        shift
        CORES="$1"
         ;;
     *)
        show_usage
        exit 1
        ;;
  esac
shift
done

    #running
    #snakemake --use-conda -s sample_test.Snakefile --configfile samp_config_test.yaml --cores 4
    
# change if statement to "if any", and/or number based

if [[ $STEP="sample" ]]; then
    paired=$(yq '.paired_end' $CONFIG);
    if [[ "${paired}" = "true" ]]; then
        snakemake --use-conda -s sample_step_paired.Snakefile --configfile $CONFIG --cores $CORES  --latency-wait 15
    else
        snakemake --use-conda -s sample_step_single.Snakefile --configfile $CONFIG --cores $CORES --latency-wait 15
    fi;
fi;

if [[ $STEP="repseqs" ]]; then
    if [[ "${asv_method}" = "dada2pe" ]]; then
        snakemake --use-conda -s repseqs_step.Snakefile dada2_pe_denoise --configfile $CONFIG --cores $CORES  --latency-wait 15
    elif [[ "${asv_method}" = "dada2se" ]]; then
        snakemake --use-conda -s srepseqs_step.Snakefile dada2_se_denoise --configfile $CONFIG --cores $CORES   --latency-wait 15
    elif [[ "${asv_method}" = "deblur" ]]; then
        snakemake --use-conda -s repseqs_step.Snakefile deblur_denoise --configfile $CONFIG --cores $CORES  --latency-wait 15
    else
        echo "ASV method not recognized"
    fi;
fi;


# #tourmaline -p config.yaml -samples
# args: config file, step, max-cores
#def sample{
#checks
#check if metadata is there

