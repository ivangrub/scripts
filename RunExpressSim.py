#!/usr/bin/env python

import os

bins = [100,500,1000,5000]

for i in bins:
	os.system('BinnedFasta.py -g mm9_norandom -b %s -l 50' % i)

for i in bins:
	os.system('bzip2 -dc CC_Rad23b.fastq.bz2 | bowtie -c -q -v 2 -k 100 -S mm9_norandom - | BinMapping.py -r - -g mm9_norandom_%s' % i)
