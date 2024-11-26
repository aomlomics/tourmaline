#!/bin/bash

# Usage: 
# scripts/initialize_dir_from_existing_tourmaline_dir.sh /path/to/old/tourmaline
#
# From the main directory of a newly cloned tourmaline directory, run this script
# to remove the test files, copy config.yaml from an existing tourmaline directory,
# then copy the data files and symlinks from the existing tourmaline directory.
# Files are only attempted to be copied if they exist.

echo "Removed:"
/bin/rm -rv 00-data/*
echo "Copied:"
[[ -f $1/config.yaml ]] && cp -av $1/config.yaml config.yaml
[[ -f $1/00-data/manifest_pe.csv ]] && cp -av $1/00-data/manifest_pe.csv 00-data/manifest_pe.csv
[[ -f $1/00-data/manifest_se.csv ]] && cp -av $1/00-data/manifest_se.csv 00-data/manifest_se.csv
[[ -f $1/00-data/metadata.tsv ]] && cp -av $1/00-data/metadata.tsv 00-data/metadata.tsv
[[ -f $1/00-data/repseqs_to_filter_dada2-pe.tsv ]] && cp -av $1/00-data/repseqs_to_filter_dada2-pe.tsv 00-data/repseqs_to_filter_dada2-pe.tsv
[[ -f $1/00-data/repseqs_to_filter_dada2-se.tsv ]] && cp -av $1/00-data/repseqs_to_filter_dada2-se.tsv 00-data/repseqs_to_filter_dada2-se.tsv
[[ -f $1/00-data/repseqs_to_filter_deblur-se.tsv ]] && cp -av $1/00-data/repseqs_to_filter_deblur-se.tsv 00-data/repseqs_to_filter_deblur-se.tsv
[[ -f $1/01-imported/refseqs.qza ]] && cp -av $1/01-imported/refseqs.qza 01-imported/refseqs.qza
[[ -f $1/01-imported/reftax.qza ]] && cp -av $1/01-imported/reftax.qza 01-imported/reftax.qza
[[ -f $1/01-imported/classifier.qza ]] && cp -av $1/01-imported/classifier.qza 01-imported/classifier.qza
