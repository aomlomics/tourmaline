#!/bin/bash
##
# @Description: Runs the snakemake pipeline Tourmaline 2.0, either all at once or modularly.
#command: ./tourmaline.sh -s sample --configfile config-01-sample.yaml --cores 4
##

function show_usage (){
    printf "Usage: $0 [options [parameters]]\n"
    printf "\n"
    printf "Options:\n"
    printf " -s|--step, sample, repseqs\n"
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

echo $STEP

case $STEP in 

    sample)
        paired=$(yq -r '.paired_end' $CONFIG);
        trim=$(yq -r '.to_trim' $CONFIG);
        echo -n "Running sample step"
        if [[ "${paired}" = true ]]; then
            if [[ "${trim}" = true ]]; then
            #snakemake --use-conda -s sample_step_paired.Snakefile --configfile $CONFIG --cores $CORES  --latency-wait 15
                snakemake --use-conda -s sample_step_paired.Snakefile trim_pe_all --configfile $CONFIG --cores $CORES  --latency-wait 15
            else
                snakemake --use-conda -s sample_step_paired.Snakefile no_trim_pe_all --configfile $CONFIG --cores $CORES  --latency-wait 15
            fi;
        else
            #snakemake --use-conda -s sample_step_single.Snakefile --configfile $CONFIG --cores $CORES --latency-wait 15
            if [[ "${trim}" = true ]]; then
                snakemake --use-conda -s sample_step_single.Snakefile trim_se_all --configfile $CONFIG --cores $CORES  --latency-wait 15
            else
                snakemake --use-conda -s sample_step_single.Snakefile no_trim_se_all --configfile $CONFIG --cores $CORES  --latency-wait 15
            fi;
        fi;
        ;;
    repseqs)
        asv_method=$(yq -r '.asv_method' $CONFIG);
        echo -n "Running repseqs step\n"
        if [[ ${asv_method} = dada2pe ]]; then
            snakemake --use-conda -s repseqs_step.Snakefile run_dada2_pe_denoise --configfile $CONFIG --cores $CORES  --latency-wait 15
        elif [[ "${asv_method}" = "dada2se" ]]; then
            snakemake --use-conda -s repseqs_step.Snakefile run_dada2_se_denoise --configfile $CONFIG --cores $CORES   --latency-wait 15
        elif [[ "${asv_method}" = "deblur" ]]; then
            snakemake --use-conda -s repseqs_step.Snakefile run_deblur_se_denoise --configfile $CONFIG --cores $CORES  --latency-wait 15
        else
            echo "ASV method not recognized"
        fi;
        ;;
    *)
        echo -n "unknown step"
        ;;
esac


# if [[ "$STEP"=sample ]]; then
#     paired=$(yq '.paired_end' $CONFIG);
#     if [[ "${paired}" = "true" ]]; then
#         snakemake --use-conda -s sample_step_paired.Snakefile --configfile $CONFIG --cores $CORES  --latency-wait 15
#     else
#         snakemake --use-conda -s sample_step_single.Snakefile --configfile $CONFIG --cores $CORES --latency-wait 15
#     fi;
# elif [[ "$STEP"="repseqs" ]]; then
#     asv_method=$(yq '.asv_method' $CONFIG);
#     if [[ "${asv_method}" = "dada2pe" ]]; then
#         snakemake --use-conda -s repseqs_step.Snakefile run_dada2_pe_denoise --configfile $CONFIG --cores $CORES  --latency-wait 15
#     elif [[ "${asv_method}" = "dada2se" ]]; then
#         snakemake --use-conda -s repseqs_step.Snakefile run_dada2_se_denoise --configfile $CONFIG --cores $CORES   --latency-wait 15
#     elif [[ "${asv_method}" = "deblur" ]]; then
#         snakemake --use-conda -s repseqs_step.Snakefile run_deblur_se_denoise --configfile $CONFIG --cores $CORES  --latency-wait 15
#     else
#         echo "ASV method not recognized"
#     fi;
# fi;



# #tourmaline -p config.yaml -samples
# args: config file, step, max-cores
#def sample{
#checks
#check if metadata is there

