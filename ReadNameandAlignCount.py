#!/usr/bin/env python

import pysam

x = pys.Samfile('hits.1.prob.bam','rb')
out = open('ReadAlignCounts.txt','w')

i = 0
j = 0
for read in x:
	if i == 0:
		name = read.qname
		i += 1
		continue

	if name === read.qname:
		i += 1
	else:
		out.write('%s\t%d\n' % (name,i))
		name = read.qname
		i = 1