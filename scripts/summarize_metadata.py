#!/usr/bin/env python

import pandas as pd
from tabulate import tabulate
import sys

# usage
usage = '''
summarize_metadata.py metadata output1_md output2_txt

'''

if len(sys.argv) < 2:
    print(usage)
    sys.exit()

# input paths
metadata = sys.argv[1]

# output paths
output1 = sys.argv[2] # 'manifest_pe.csv'
output2 = sys.argv[3] # 'manifest_se.csv'

df = pd.read_csv(metadata, sep='\t')
cols = df.columns
df2 = pd.DataFrame(columns =[0,1], index=cols)
for col in cols:
    if col in df.columns:
        vc = df[col].value_counts()
        if vc.index.shape == (0,):
            df2.loc[col, 0] = '(no values in column)'
            df2.loc[col, 1] = '--'
        else:
            df2.loc[col, 0] = vc.index[0]
            df2.loc[col, 1] = vc.values[0]
    else:
        df2.loc[col, 0] = '(column not provided)'
        df2.loc[col, 1] = '--'
df2.columns = ['Most common value', 'Count']
df2.index.name = 'Column name'
outstr = tabulate(df2, tablefmt="pipe", headers="keys")
with open(output1, 'w') as target:
    target.write(outstr)
    target.write('\n')
with open(output2, 'w') as target:
    [target.write('%s\n' % i) for i in cols]
