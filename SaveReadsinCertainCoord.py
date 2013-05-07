#!/usr/bin/env python

import pysam as pys
import numpy as np
import argparse
import sys

def GetReads(chr,coord,lcoord,rcoord,sam,out):
	out.write('%s\n' % args.r)
	i = 0
	for read in sam:
		if ref[read.tid] == chr and read.pos >= lcoord and read.pos <= rcoord:
			try:
				frac_str = read.opt('XP')
			except KeyError:
				frac_str = '1'
			out.write('%s\t%s:%d-%d\t%s\n' % (read.qname,chr,read.pos,read.pos+50,frac_str))
		if (i % 1000000) == 0:
			print '%s reads' % i
		i += 1

	out.close()
	sam.close()

parser = argparse.ArgumentParser(description='Return reads that map to a coordinate range')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = None)
parser.add_argument('-c',help='The coordinates in chrX:1-10 format',default =None)
parser.add_argument('-p',help='Input a peak file',default = None)
parser.add_argument('-o',help='Output tag that will be assigned to MultiplyAligned file',default='output')
args = parser.parse_args()

if (args.c == None or args.p == None) or args.r == None:
	sys.exit('Enter the coordinates or peak file and the BAM file')

if (args.c != None):	
	s = args.c.split(':')
	chr = s[0]
	coord = s[1].split('-')
	lcoord = int(coord[0])
	rcoord = int(coord[1])
	coordinates = '%s_%d_%d' % (chr,lcoord,rcoord)	
	outfile =open('%s_MultiplyAligned_%s.txt' % (args.o,coordinates),'w')
	samfile = pys.Samfile('%s' % args.r,'rb')
	ref = samfile.references	
	GetReads(chr,coord,lcoord,rcoord,samfile,outfile)
else:
	peaks = open('%s' % args.p,'r')
	for line in peaks:
		s = line.strip().split('\t')
		coordinates = '_'.join(s[0:2])
		chr = s[0]
		lcoord = int(s[1])
		rcoord = int(s[2])
		outfile =open('%s_MultiplyAligned_%s.txt' % (args.o,coordinates),'w')
		samfile = pys.Samfile('%s' % args.r,'rb')
		ref = samfile.references	
		GetReads(chr,coord,lcoord,rcoord,samfile,outfile)
