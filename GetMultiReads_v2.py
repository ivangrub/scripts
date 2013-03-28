#!/usr/bin/env python

import pysam as pys

samfile = pys.Samfile('BG_dm3.allalign.sam','r')
outfile = pys.Samfile('BG_dm3_multireads.sam','w',template = samfile)
i = 0
j = 0
for read in samfile:
	if i > 0: 
		if x.seq == read.seq:
			outfile.write(x)
			outfile.write(read)
			j += 1
	x = read
	i += 1
	if i % 1e6 == 0:
		print '%d reads processed' % i

print 'Number of multiple alignments: %d' % j		
outfile.close()
samfile.close()