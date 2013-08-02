#!/usr/bin/env python

import pysam as pys
import numpy as np
import argparse
#import yappi

def GetLikelihood(seq,ind):
	frac_str = seq.opt('XP')
	frac = float(frac_str)
	hist[ind][int(frac*100)-1] += 1

def FindOverlap(seq,coords,chromo,ref):
	chrom = ref[seq.tid]
	st = seq.pos
	a = chrom == chromo
	c = st >= coords[:,0]
	d = st <= coords[:,1]
	cboth = c == d
	out = np.where(np.logical_and(a,cboth))[0]
	
	return out
	

parser = argparse.ArgumentParser(description='Project SAM/BAM format reads back onto the appropriate genome and call significant peaks.')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = '')
parser.add_argument('-p',help = 'The peak file in BED format',default = '')
parser.add_argument('-o',help = 'Prefix of output file containing the histogram',default = 'output')
args = parser.parse_args()

bam = pys.Samfile(args.r,'rb')
peak = open(args.p,'r')
hist_out = open(args.o,'w')

refdict = bam.references
coord = np.array([])
chr = np.array([[]])
name = []
for line in peak:
	s = line.strip().split('\t')
	name.append(":".join(s))
	chr = np.append(chr,[[s[0]]])
	coord = np.append(coord,[float(s[1]),float(s[2])])

peak.close()
x = len(coord)
crd = np.reshape(coord,(x/2,2))
hist = np.zeros((x/2,100))
#yappi.start()
i = 0
for read in bam:
	i += 1
	if (i % 1000000 == 0):
		print 'Processed %d Alignments' % i

	index = FindOverlap(read,crd,chr,refdict)
	if len(index) > 0:
		GetLikelihood(read,index[0])
		
	#yappi.print_stats()
	


print 'Printing Text file'
for i in xrange(np.shape(hist)[0]):
	hist_out.write('%s\t%s\n' % (name[i],"\t".join(map(str,hist[i]))))

hist_out.close()

bam.close()
