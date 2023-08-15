#!/usr/bin/env python

import pandas as pd
import sys
from qiime2 import Artifact

# usage
usage = '''
filter_taxonomy.py taxonomy repseqs taxonomytsv taxonomyqza
'''

if len(sys.argv) < 3:
    print(usage)
    sys.exit()

# input paths
taxonomy = sys.argv[1]
repseqs = sys.argv[2]

# output paths
taxonomytsv = sys.argv[3] #
taxonomyqza = sys.argv[4]

df_taxonomy = pd.read_csv(taxonomy, index_col=0, sep='\t')
df_repseqs = pd.read_csv(repseqs, header=None, index_col=0, sep='\t')
keep_ids = df_repseqs.index
df_taxonomy_filtered = df_taxonomy.loc[list(keep_ids)]
df_taxonomy_filtered.to_csv(taxonomytsv, sep='\t')
artifact_taxonomy_filtered = Artifact.import_data('FeatureData[Taxonomy]', df_taxonomy_filtered)
artifact_taxonomy_filtered.save(taxonomyqza)