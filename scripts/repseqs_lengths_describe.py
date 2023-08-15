#!/usr/bin/env python

import pandas as pd
import sys
from tabulate import tabulate

# usage
usage = '''
repseqs_lengths_describe.py repseqs_lengths_tsv repseqs_lengths_md

'''

if len(sys.argv) < 1:
    print(usage)
    sys.exit()

# input paths
IN = sys.argv[1]

# output paths
OUT = sys.argv[2] # 'repseqs_lengths_describe.md'

s = pd.read_csv(IN, header=None, index_col=0, sep='\t')
t = s.describe()
outstr = tabulate(t.iloc[1:], tablefmt="pipe", headers=['Statistic (n=%s)' % t.iloc[0].values[0].astype(int), 'Sequence length'])
with open(OUT, 'w') as target:
    target.write(outstr)
    target.write('\n')

