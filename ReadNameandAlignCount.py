#!/usr/bin/env python

import pysam as pys

x = pys.Samfile('hits.1.prob.bam','rb')
out = open('ReadAlignCounts.txt','w')

i = 0
j = 0
for read in x:
	if i == 0:
		name = read.qname
		i = 1
		j += 1
		continue

	if name == read.qname:
		j += 1
	else:
		out.write('%s\t%d\n' % (name,j))
		name = read.qname
		j = 1

out.write('%s\t%d\n' % (name,j))
x.close()
out.close()