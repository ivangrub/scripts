#!/usr/bin/env python

import pysam as pys

samfile = pys.Samfile('PI-k100/Post-eXpress-PI-k100.dm3.bam','rb')
fasta = open('PI-k100/Multi-seq-Left.fa','w')

ref = samfile.references
chrom = 'chr2L'
coord = [21473000,21473500]
k = 0
i = 0
for read in samfile:
	if ref[read.tid] == chrom and read.pos > coord[0] and read.pos < coord[1]:
		fasta.write('>%s\n%s\n' % (read.qname,read.seq))

	k += 1
	if k % 1e6 == 0:
		print 'read %d reads' % k
		
fasta.close()