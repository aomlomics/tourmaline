#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
from tabulate import tabulate
from qiime2 import Artifact
import seaborn as sns
import matplotlib.pyplot as plt

# usage
usage = '''
plot_repseq_properties.py repseqs_lengths_tsv aligned_repseqs_gaps aligned_repseqs_outliers taxonomy_qza table_qza repseqs_properties repseqs_properties_describe repseqs_properties_pdf outliers_tsv

'''

if len(sys.argv) < 8:
    print(usage)
    sys.exit()

# input paths
lengths = sys.argv[1]
gaps = sys.argv[2]
outliers = sys.argv[3]
taxonomy = sys.argv[4]
table = sys.argv[5]

# output paths
proptsv = sys.argv[6] #
propdescribe = sys.argv[7]
proppdf = sys.argv[8]
outliersforqza = sys.argv[9] 

lengths = pd.read_csv(lengths, names=['length'], index_col=0, sep='\t')
gaps = pd.read_csv(gaps, names=['gaps'], index_col=0, sep='\t')
outliers = pd.read_csv(outliers, names=['outlier'], index_col=0, sep='\t')
taxonomy = Artifact.load(taxonomy)
taxonomydf = taxonomy.view(view_type=pd.DataFrame)
taxonomydf['level_1'] = [x.split(';')[0] for x in taxonomydf['Taxon']]
table = Artifact.load(table)
tabledf = table.view(view_type=pd.DataFrame)
merged = pd.merge(lengths, gaps, left_index=True, right_index=True, how='outer')
merged = pd.merge(merged, outliers, left_index=True, right_index=True, how='outer')
merged = pd.merge(merged, taxonomydf['Taxon'], left_index=True, right_index=True, how='outer')
merged = pd.merge(merged, taxonomydf['level_1'], left_index=True, right_index=True, how='outer')
merged = pd.merge(merged, tabledf.sum().rename('observations'), left_index=True, right_index=True, how='outer')
merged.columns = ['length', 'gaps', 'outlier', 'taxonomy', 'taxonomy_level_1', 'observations']
merged.index.name = 'featureid'
merged['log10(observations)'] = [np.log10(x) for x in merged['observations']]
merged.sort_values('log10(observations)', ascending=False, inplace=True)
merged.to_csv(proptsv, index=True, sep='\t')
t = merged.describe()
tcolumns = t.columns
tcolumns = tcolumns.insert(0, 'Statistic (n=%s)' % t.iloc[0].values[0].astype(int))
outstr = tabulate(t.iloc[1:], tablefmt="pipe", headers=tcolumns)
with open(propdescribe, 'w') as target:
    target.write(outstr)
    target.write('\n')
g = sns.relplot(data=merged, x='length', y='gaps', col='outlier', hue='taxonomy_level_1', size='log10(observations)', sizes=(1,500), edgecolor = 'none', alpha=0.7)
g.set_axis_labels('length (bp) not including gaps', 'gaps (bp) in multiple sequence alignment')
plt.savefig(proppdf, bbox_inches='tight')
outliers.columns = ['Outlier']
outliers.index.name = 'Feature ID'
outliers = outliers*1
outliers.to_csv(outliersforqza, index=True, sep='\t')

