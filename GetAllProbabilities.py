#!/usr/bin/env python

import pysam as pys
import argparse

parser = argparse.ArgumentParser(description='Get all of the posterior probabilities of read alignments')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = '')
parser.add_argument('-o',help ='Name of the output file',default='output.txt')
args = parser.parse_args()

alignprob = open(args.o,'w')
align = pys.Samfile(args.r,'rb')

for read in align:
	try:
		frac_str = read.opt('XP')
		alignprob.write('%s\t%s\n' %(read.qname,frac_str))
	except KeyError:
		continue

alignprob.close()
align.close()