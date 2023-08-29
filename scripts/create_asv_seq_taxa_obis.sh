#!/bin/bash
##
# @Description: Takes Qiime2 generated feature table, representative sequences, and taxonomy Artifact files to generate a file formatted for OBIS conversion.

##

function show_usage (){
    printf "Usage: $0 [options [parameters]]\n"
    printf "\n"
    printf "Options:\n"
    printf " -f|--feature-table, Feature table .qza\n"
    printf " -t|--taxonomy, taxonomy.qza\n"
    printf " -r|--repseqs, repseqs.qza\n"
    printf " -o|--output, output tsv\n"
    printf " -h|--help, Print help\n"

return 0
}


while [ ! -z "$1" ]; do
  case "$1" in
     --feature-table|-f)
         shift
         FEAT="$1"
         ;;
     --taxonomy|-t)
         shift
         TAX="$1"
         ;;
     --repseqs|-r)
        shift
        REP="$1"
         ;;
    --output|-o)
        shift
        OUT="$1"
         ;;
     *)
        show_usage
        exit 1
        ;;
  esac
shift
done

qiime feature-table transpose \
  --i-table $FEAT \
  --o-transposed-feature-table transposed-table.qza

qiime metadata tabulate \
  --m-input-file $REP \
  --m-input-file $TAX \
  --m-input-file transposed-table.qza \
  --o-visualization merged-data.qzv

qiime tools export \
  --input-path merged-data.qzv \
  --output-path $OUT

mv $OUT/metadata.tsv temp 
rm -r $OUT
sed -e '2d' temp | sed "1 s|id\t|featureid\t|" | sed "1 s|Taxon|taxonomy|" | sed "1 s|Sequence|sequence|" > $OUT

rm temp transposed-table.qza merged-data.qzv