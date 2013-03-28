#!/usr/bin/env python

import pysam as pys
import numpy as np

first = pys.Samfile('Sim_TBP.900bp.k100_converted.sorted.bam','rb')
second = pys.Samfile('sim_900bp_wB/hits.1.prob.bam','rb')

f_id = {}
f_pos = {}
print 'Loading Dictionary'
for read in first:
	try:
		f_id[read.qname].append(read.tid)
		f_pos[read.qname].append(read.pos)
	except KeyError:
		f_id[read.qname] = [read.tid]
		f_pos[read.qname] = [read.pos]


print 'Checking express processed reads'
ref = second.references
for read in second:
	s = read.qname

	try:
		fpos = np.array(f_pos[s])
		fid = np.array(f_id[s])
		lpos = fpos == read.pos
		lid = fid == read.tid
		if np.logical_xor(lpos,lid).any() == False:
			print read.pos,ref[read.tid],read.tid	
	except KeyError:
		print read.pos,ref[read.tid],read.tid

first.close()
second.close()
