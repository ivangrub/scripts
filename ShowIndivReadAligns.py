#!/usr/bin/env python

import pysam as pys

bam = pys.Samfile('','rb')
out = open('alignof_16604.txt','w')
name = '16604'

for read in bam:
	if name == read.qname:
		out.write('%s\t%s' %(read.qname,read.pos))
		
