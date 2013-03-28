#!/usr/bin/env python

import pysam as pys
import numpy as np

region = 'chr2L:21,415,513-21,545,302'
file = 'Fly_Pol2.k500.dm3.bam'
bam = pys.Samfile('%s' % file,'rb')
out = open('UniqueinRegion_%s_%s.txt' %(region,file[:-4]),'w')

reg = region.split(':')
chr = reg[0]
nocomma = reg[1].replace(',','')
coord = reg.split('-')
ref = bam.references
i = 0
count = 0
chrom = []
pos = []
for read in bam:
	if i == 0:
		name = read.qname
		count += 1
		pos.append(read.pos)
		chrom.append(ref[read.tid])
		i += 1
		continue
	if name == read.qname:
		count += 1
		pos.append(read.pos)
		chrom.append(ref[read.tid])
		continue
	else:
		ch = np.array(chrom)
		st = np.array(pos)
		if np.logical_and(np.logical_and(chr,ch),np.logical_and(st >= coord[0],st<= coord[1])) and count == 1:
			out.write('%s\n' % (name))
		count = 1
		name = read.qname
		chrom = []
		pos = []

out.close()
bam.close()