#!/usr/bin/env python

import pandas as pd
import sys
from tabulate import tabulate

# usage
usage = '''
describe_fastq_counts.py fastq_counts_tsv output1_md

'''

if len(sys.argv) < 1:
    print(usage)
    sys.exit()

# input paths
IN = sys.argv[1]

# output paths
OUT = sys.argv[2] # 'fastq_counts_describe.md'


s = pd.read_csv(IN, index_col=0, sep='\t')
t = s.describe()
outstr = tabulate(pd.DataFrame(t.iloc[1:,0]), tablefmt="pipe", headers=['Statistic (n=%s)' % t.iloc[0,0].astype(int), 'Fastq sequences per sample'])
with open(OUT, 'w') as target:
    target.write(outstr)
    target.write('\n')

