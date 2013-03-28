#!/usr/bin/env python

import pysam as pys
import numpy as np

region = 'chr2L:21,415,513-21,545,302'
file = 'Fly_Pol2.m1.dm3.bam'
bam = pys.Samfile('%s' % file,'rb')
out = open('UniqueinRegion_%s_%s.txt' %(region,file[:-4]),'w')

reg = region.split(':')
chr = reg[0]
nocomma = reg[1].replace(',','')
coord = nocomma.split('-')
ref = bam.references
i = 0
count = 0
chrom = []
pos = []
for read in bam:
	if chr == ref[read.tid] and read.pos >= coord[0] and read.pos <= coord[1]:
		out.write('%s\n' % read.qname)
	i += 1
	if i % 1000000 == 0:
		print '%d reads processed' % i

out.close()
bam.close()