#!/usr/bin/env python

import pysam as pys

dm3 = pys.Samfile('Pol2.dm3.bam','rb')
dm3out = pys.Samfile('Pol2.dm3.unmapped.sam','w',template = dm3)

i = 0
for read in dm3:
	if read.is_unmapped:
		dm3out.write(read)
	i += 1
	if i % 1e6 == 0:	
		print '%d reads processed' % i

dm3out.close()
dm3.close()