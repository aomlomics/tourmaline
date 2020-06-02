#!/bin/bash

# Usage: 
# scripts/copy_data_files_to_new_tourmaline_dir.sh /path/to/old/tourmaline
#
# From the main directory of a newly cloned tourmaline directory, run this
# script to first remove the test files, then copy the data files and symlinks
# of an existing tourmaline directory on your system.

cp -a $1/config.yaml .
rsync -a $1/Snakefile .
/bin/rm 00-data/*
cp -a $1/00-data/manifest_pe.csv 00-data/manifest_pe.csv
cp -a $1/00-data/manifest_se.csv 00-data/manifest_se.csv
cp -a $1/00-data/metadata.tsv 00-data/metadata.tsv
cp -a $1/00-data/refseqs.fna 00-data/refseqs.fna
cp -a $1/00-data/reftax.tsv 00-data/reftax.tsv
mkdir 01-imported
cp -a $1/01-imported/refseqs.qza 01-imported/refseqs.qza
cp -a $1/01-imported/reftax.qza 01-imported/reftax.qza
