#!/usr/bin/env python

import pysam as pys

out = open('ReadProbability.txt','w')
reads = pys.Samfile('hits.1.prob.bam','rb')

for r in reads:
	prob = r.opt('XP')
	out.write('%s\t%s\n' % (r.qname,prob))