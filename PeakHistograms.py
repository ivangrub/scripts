#!/usr/bin/env python

import pysam as pys
import numpy as np
import argparse


parser = argparse.ArgumentParser(description='Project SAM/BAM format reads back onto the appropriate genome and call significant peaks.')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = '')
parser.add_argument('-p',help = 'The peak file in BED format',default = '')
parser.add_argument('-o',help = 'Prefix of output file containing the histogram',default = 'output')
args = parser.parse_args()

bam = pys.Samfile(args.r,'rb')
peak = open(args.p,'r')
hist_out = open(args.o,'w')

coord = np.array([])
chr = np.array([[]])
name = []
for line in peak:
	s = line.strip().split('\t')
	name.append(":".join(s))
	chr = np.append(chr,[[s[0]]])
	coord = np.append(coord,[float(s[1]),float(s[2])])

x = len(coord)
crd = np.reshape(coord,(x/2,2))
hist = np.zeros((x/2,100))

i = 0
for read in bam:
	i += 1

	chrom = bam.getrname(read.tid)
	st = read.pos
	a = chrom == chr
	c = st >= crd[:,0]
	d = st <= crd[:,1]
	cboth = np.logical_and(c,d)

	if (i % 1000000) == 0:
		print "%d reads processed" % i
	try:
		index = np.where(np.logical_and(a,cboth))[0][0]
	except IndexError:
		continue
	
	frac_str = read.opt('XP')
	frac = float(frac_str)
	try:
		likeli = int("%.2f" % round(frac))
		hist[index][int(likeli*100)] += 1
	except ValueError:
		hist[index][99] += 1


print 'Printing Text file'
for i in xrange(np.shape(hist)[0]):
	hist_out.write('%s\t%s\n' % (name[i],"\t".join(map(str,hist[i]))))

hist_out.close()
peak.close()
bam.close()
