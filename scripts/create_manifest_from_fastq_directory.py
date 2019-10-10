#!/usr/bin/env python

import sys
import os
import glob
import re

# usage
usage = '''
create_manifest_from_fastq_directory.py FASTQ_DIR MANIFEST_PE MANIFEST_SE

    FASTQ_DIR - full path of directory containing fastq.gz files
        of the form CN18SESPkoa_SC36_S80_L001_R1_001.fastq.
	MANIFEST_PE - output path of manifest_pe.csv
	MANIFEST_SE - output path of manifest_se.csv

    This script makes the following assumptions:
        - First bit after sample name is "_S[0-9]{1,3}""
        - Only Lane 1 data is present "_L001"
        - R1 and R2 files are both present "_R1"/_R2"
        - Only one file "_001" for each of R1 and R2
'''

if len(sys.argv) < 2:
    print(usage)
    sys.exit()

# input paths
path_fastq = sys.argv[1]

# output paths
path_manifest_pe = sys.argv[2] # 'manifest_pe.csv'
path_manifest_se = sys.argv[3] # 'manifest_se.csv'

# list of full fastq.gz paths
list_fastq_full = sorted(glob.glob(os.path.join(path_fastq, '*.fastq.gz')))

# list of sample names (one per path)
list_fastq = [x.split('/')[-1] for x in list_fastq_full]
list_samples = [re.sub('_S[0-9]{1,3}_L001_R[12]_001.fastq.gz', '', x) for x in list_fastq]

# list of forward and reverse (one of each per sample)
list_fwdrev = ['forward', 'reverse'] * int(len(list_samples)/2)

# write manifest_pe
with open(path_manifest_pe, 'w') as target:
	target.write('sample-id,absolute-filepath,direction\n')
	for i in range(len(list_samples)):
		target.write('%s,%s,%s\n' % (list_samples[i], list_fastq_full[i], list_fwdrev[i]))

# get alternate rows for manifest_se
list_samples_se = list_samples[::2]
list_fastq_full_se = list_fastq_full[::2]
list_fwdrev_se = list_fwdrev[::2]

# write manifest_se
with open(path_manifest_se, 'w') as target:
	target.write('sample-id,absolute-filepath,direction\n')
	for i in range(len(list_samples_se)):
		target.write('%s,%s,%s\n' % (list_samples_se[i], list_fastq_full_se[i], list_fwdrev_se[i]))
