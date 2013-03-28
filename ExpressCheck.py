#!/usr/bin/env python

import os

path = os.environ['EXPRESS_FILES']

bins = [100,500,1000,5000]

for i in bins:
	os.system('BinnedFasta.py -g mm9_norandom.fa -l 50 -b %d' % i)
	os.system('BinMapping.py -r TY_ES_TBP.mm9_norandom.k100.bam -b %d -o TBP_bin%d' % (i,i))

neigh = [1000,5000,10000]
for i in bins:
	for j in neigh:
		n = 0.5*j/i
		os.system('express -o express_out_bin%d --num-neighbors %d -B 1 %s/eXpress_%dbp_50.mm9_norandom.fa TBP_bin%d_%dneigh_converted.bam' % (i,n,path,i,i,n))
		os.system('eXpress2wiggle.py -r TBP_bin%d_converted.bam -g mm9_norandom -o PostExpress_TBP_bin%d_%dneigh' % (i,i,n))