#!/bin/bash
##
# @Description: Runs the snakemake pipeline Tourmaline 2.0, either all at once or modularly.
#command: ./tourmaline.sh -s sample --configfile config-01-sample.yaml --cores 4
##

function show_usage (){
    printf "Usage: $0 [options [parameters]]\n"
    printf "\n"
    printf "Options:\n"
    printf " -s|--step, samples\n"
    printf " -c|--configfile, taxonomy.qza\n"
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
    
if [[ $STEP="samples" ]]; then
    run_name=$(yq '.run_name' $CONFIG);
    paired=$(yq '.paired_end' $CONFIG);
    if [[ "${paired}" = "true" ]]; then
        snakemake --use-conda -s sample_test_paired.Snakefile --configfile $CONFIG --cores $CORES  --latency-wait 15
    else
        snakemake --use-conda -s sample_test_single.Snakefile --configfile $CONFIG --cores $CORES --conda-frontend conda
    fi;
fi;


# #tourmaline -p config.yaml -samples
# args: config file, step, max-cores
#def sample{
#checks
#check if metadata is there

