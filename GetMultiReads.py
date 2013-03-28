#!/usr/bin/env python

import pysam as pys

samfile = pys.Samfile('GRO-seq/Post-eXpress-GRO-all.dm3.bam','rb')
outfile = pys.Samfile('GRO-seq/Multireads.sam','w',template = samfile)

i = 0
j = 0
k = 0
for read in samfile:
	if i > 0: 
		if x.qname == read.qname:
			if k == 0:
				outfile.write(x)
				outfile.write(read)
				k = 1
			else:
				outfile.write(read)
			j += 1
		else:
			k = 0
	x = read
	i += 1
	if i % 1e6 == 0:
		print '%d reads processed' % i

print 'Number of multiple alignments: %d' % j		
outfile.close()
samfile.close()