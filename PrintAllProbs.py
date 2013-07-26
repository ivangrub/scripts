#!/usr/bin/env python

import pysam as pys

out = open('ReadProbability.txt','w')
reads = pys.Samfile('CSEM_Fly_Pol2.k500.sim_hisrepeat.bam','rb')

for r in reads:
	try:
		prob = r.opt('ZW')
	except KeyError:
		continue
	out.write('%s\t%s\n' % (r.qname,prob))