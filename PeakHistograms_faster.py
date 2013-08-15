#!/usr/bin/env python

import pysam as pys
import numpy as np
import argparse
import cProfile, pstats, io
#import yappi

def GetLikelihood(seq,ind):
	frac_str = seq.opt('XP')
	frac = float(frac_str)
	index = int(frac*100)-1
	try:
		hist[ind][index] += 1
	except IndexError:
		hist[ind][index] += 1

def FindOverlap(seq,peaks,ref):
	pr.enable()
	chrom = ref[seq.tid]
	st = seq.pos
	for val in peaks.values():
		if (val[1] == chrom and st >= val[2] and st <= val[3]):
			return val[0]
	return False
	


parser = argparse.ArgumentParser(description='Project SAM/BAM format reads back onto the appropriate genome and call significant peaks.')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = '')
parser.add_argument('-p',help = 'The peak file in BED format',default = '')
parser.add_argument('-o',help = 'Prefix of output file containing the histogram',default = 'output')
args = parser.parse_args()

bam = pys.Samfile(args.r,'rb')
peak = open(args.p,'r')
hist_out = open(args.o,'w')

refdict = bam.references
peakdict = {}
#coord = np.array([])
#chr = np.array([[]])
name = []
#counter = np.array([])
i = 0
for line in peak:
	s = line.strip().split('\t')
	name.append(":".join(s))
	peakdict[name[i]] = [i,s[0],float(s[1]),float(s[2])]
	#chr = np.append(chr,[[s[0]]])
	#coord = np.append(coord,[float(s[1]),float(s[2])])
	#counter = np.append(counter,[i])
	i += 1

peak.close()
x = len(peakdict.keys())
#crd = np.reshape(coord,(x/2,2))
hist = np.zeros((x,100))

#yappi.start()
i = 0
pr = cProfile.Profile()

for read in bam:

	i += 1
	if (i % 1000000 == 0):
		print 'Processed %d Alignments' % i

	index = FindOverlap(read,peakdict,refdict)

	if index != False:
		GetLikelihood(read,index)
		
	#yappi.print_stats()
	
	


print 'Printing Text file'
for i in xrange(np.shape(hist)[0]):
	hist_out.write('%s\t%s\n' % (name[i],"\t".join(map(str,hist[i]))))

hist_out.close()

bam.close()
