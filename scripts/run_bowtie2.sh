#!/bin/bash

## Run bowtie2 on repseqs.fasta
### this script is run as follows
# bowtie2.sh -o <out_dir_for_anacapa_QC_run> -d <database_directory_index> -t temp_dir -r <repseqs.fasta>

OUT=""
T_DIR=""
DB=""
REPSEQS=""

while getopts "o:t:d:r:" opt; do
    case $opt in
        o) OUT="$OPTARG" # path to desired bowtie2 output
        ;;
        t) T_DIR="$OPTARG" # path to temp output directory
        ;;
        d) DB="$OPTARG"  # path to reference database bowtie2 index
        ;;
        r) REPSEQS="$OPTARG"  # path to repseqs.fasta
        ;;
    esac
done


######################################
# Submit single read and un merged paired read ASV's to Bowtie2
######################################
mkdir -p ${T_DIR}

echo "Run Bowtie2 on ASV fasta file"
# find best taxonomic hits for single reads in using bowtie2's global alignment mode (end-to-end) first.  It only considers full length alighnments, and recovers up to -k returns (default = 100). The reads not assigned in global mode are then run in local alighnment mode. -p is the number of threads --no-hd and --no-sq suppress the header lines --no-unal means do not add unaligned reads to the sam file (output). --very-sensitive is a preset option in bowtie2 is slow but designed to be more accurate and sensitive

if [[ -e "${REPSEQS}" ]];
  then
    echo ""
    echo "${REPSEQS} exists"
    echo "end-to-end"
    bowtie2 -x ${DB}  -f -U ${REPSEQS} -S ${T_DIR}/end_to_end.sam --no-hd --no-sq --very-sensitive --end-to-end --no-unal -p 120 -k 100 --un ${T_DIR}/end_to_end_reject.fasta
    date
    #single read local mode
    # reads that do not get hits in global alignment mode are run through local alignment.  This mode does not require the hit to be the entire length of the query.  This includes partial matches to full length references, or full matches to short references.
    echo "local"
    bowtie2 -x ${DB}  -f -U ${T_DIR}/end_to_end_reject.fasta -S ${T_DIR}/local.sam --no-hd --no-sq --very-sensitive --local --no-unal -p 120 -k 100 --un ${T_DIR}/end_to_end_and_local_reject.fasta
    date
  else
    echo ""
    echo "${REPSEQS} does not exists"
  fi


######################################
# concatenate bowtie2 tables 
######################################

### concat all of the sam files for blca
cat ${T_DIR}/*.sam > ${OUT}
rm ${T_DIR}/end_to_end.sam
rm ${T_DIR}/local.sam