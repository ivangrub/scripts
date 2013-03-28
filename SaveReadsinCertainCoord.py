#!/usr/bin/env python

import pysam as pys
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(description='Return reads that map to a coordinate range')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = None)
parser.add_argument('-c',help='The coordinates in chrX:1-10 format',default =None)
parser.add_argument('-o',help='Output tag that will be assigned to MultiplyAligned file',default='output')
args = parser.parse_args()

if args.c == None or args.r == None:
	sys.exit('Enter the coordinates or BAM file')
	
s = args.c.split(':')
chr = s[0]
coord = s[1].split('-')
lcoord = int(coord[0])
rcoord = int(coord[1])


outfile =open('MultiplyAligned_%s.txt' % args.c,'w')
	
samfile = pys.Samfile('%s' % args.r,'rb')
ref = samfile.references
	
i = 0
outfile.write('%s\n' % args.r)
for read in samfile:
	if ref[read.tid] == chr and read.pos >= lcoord and read.pos <= rcoord:
		try:
			frac_str = read.opt('XP')
		except KeyError:
			frac_str = '1'
		outfile.write('%s\t%s:%d-%d\t%s\n' % (read.qname,chr,read.pos,read.pos+50,frac_str))
	if (i % 1000000) == 0:
		print '%s reads' % i
	i += 1

outfile.close()
samfile.close()
