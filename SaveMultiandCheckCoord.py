#!/usr/bin/env python

import pysam as pys
import numpy as np

samfile = pys.Samfile('/Users/ivang/Desktop/Input-all.bam','rb')

outfile =open('Input-all-100kb.ReadsAndLocations.txt','w')
outfile.write('ID\tChr\tPosition\tSequence\tStrand\n')

counts = open('Input-all-100kb.Counts.txt','w')
counts.write('ID\t# Chr2L\t# ~Chr2L\t# In Coordinates\n')

ref = samfile.references
chrom = 'chr2L'
coord = [[21315778,21641719]]

i = 0
k = 0
r = []
y = []
for read in samfile:
	if read.is_unmapped:
		continue
	if k == 0:
		r.append(ref[read.tid])
		y.append(read.pos)
		k = 1
	elif read.qname != x.qname:
		if chrom in r:
			if x.is_reverse:
				strand = '-'
			else:
				strand = '+'
			t = True
			ind2 = 0
			while t:
				try:
					ind = r[ind2:].index(chrom)
					ind = ind + ind2
					for z in coord:
						if y[ind] > z[0] and y[ind] < z[1]:
							for zz in range(len(y)):
								outfile.write('%s\t%s\t%s\t%s\t%s\n' % (x.qname,r[zz],y[zz],x.seq,strand))
							counts.write('%s\t%d\t%d\t%d\n' % (x.qname,r.count(chrom),len(r)-r.count(chrom),sum(np.logical_and(np.greater(y,z[0]),np.less(y,z[1])))))	
						t = False
						break
				except ValueError:
					t = False
					break
				ind2 = ind + 1
		k = 0
		r = []
		y = []
		continue
	if k == 1:
		r.append(ref[read.tid])
		y.append(read.pos)
	i+=1
	x = read
	if i % 1e6 == 0:
		print '%d locations processed' % i
	
outfile.close()
samfile.close()
counts.close()
