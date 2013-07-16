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
hist = open(args.o,'w')

coord = np.zeros((50000,2))
chr = np.empty([50000,1])
name = []
i = 0
for line in peak:
	s = line.strip().split('\t')
	name.append(":".join(s))
	chr = np.append(chr,[[s[0]]])
	coord[i] = [int(s[1]),int(s[2])]
	i += 1

chr = chr[0:i][:]
coord = coord[0:i][:]
print chr,coord
hist = np.matrix((i,100))

for read in bam:
	chrom = read.getrname
	st = read.pos
	a = chrom == chr
	c = st > coord[:][0]
	d = st < coord[:][1]
	cboth = np.logical_and(c,d)
	index = np.where(np.logical_and(a,cboth))

	frac_str = read.opt('XP')
	frac = float(frac_str)

	likeli = int("%.2f" % round(frac))
	hist[index][int(likeli*100)] += 1

for i in xrange(np.shape(hist)[0]):
	hist.write('%s\t%s\n' % (name[i],"\t".join(map(str,hist[i]))))

hist.close()
peak.close()
bam.close()
