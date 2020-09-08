#!/bin/bash

# Usage: 
# scripts/initialize_dir_from_existing_tourmaline_dir.sh /path/to/old/tourmaline
#
# From the main directory of a newly cloned tourmaline directory, run this
# script to copy to config.yaml and Snakefile from an existing tourmaline
# directory, remove the test files, then copy the data files and symlinks
# from the existing tourmaline directory.

cp -a $1/config.yaml .
rsync -a $1/Snakefile .
/bin/rm -r 00-data/*
cp -a $1/00-data/manifest_pe.csv 00-data/manifest_pe.csv
cp -a $1/00-data/manifest_se.csv 00-data/manifest_se.csv
cp -a $1/00-data/metadata.tsv 00-data/metadata.tsv
cp -a $1/01-imported/refseqs.qza 01-imported/refseqs.qza
cp -a $1/01-imported/reftax.qza 01-imported/reftax.qza
cp -a $1/01-imported/classifier.qza 01-imported/classifier.qza
