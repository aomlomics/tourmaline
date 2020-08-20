#!/usr/bin/env python

from sys import argv
import sys
import pandas as pd

# usage
usage = '''
create_manifest_from_metadata.py METADATA MANIFEST MANIFEST_PE MANIFEST_SE

    METADATA - input metadata/mapping file with sample names in first column;
        tab-delimited format (metadata.tsv)
    MANIFEST - input manifest containing all sequences files, forward and reverse;
        paths must be absolute filepaths; csv format (manifest.csv)
    MANIFEST_PE - output manifest file containing only sequence files found in 
        metadata file, forward and reverse; csv format (manifest_pe.csv)
    MANIFEST_SE - output manifest file containing only sequence files found in
        metadata file, forward only; csv format (manifest_se.csv)
'''

if len(argv) < 4:
    print(usage)
    sys.exit()

# input paths
path_metadata = argv[1] # metadata.tsv
path_manifest = argv[2] # manifest.csv

# output paths
path_manifest_pe = argv[3] # manifest_pe.csv
path_manifest_se = argv[4] # manifest_se.csv

# import manifest
df_manifest = pd.read_csv(path_manifest, index_col=0)

# import metadata
df_metadata = pd.read_csv(path_metadata, index_col=0, sep='\t')

# select manifest rows matching metadata sample ids (paired-end and single-end)
df_manifest_pe = df_manifest.loc[df_metadata.index].dropna()
df_manifest_se = df_manifest[df_manifest.direction == 'forward'].loc[df_metadata.index].dropna()

# write to csv
df_manifest_pe.to_csv(path_manifest_pe, index_label='sample-id')
df_manifest_se.to_csv(path_manifest_se, index_label='sample-id')
