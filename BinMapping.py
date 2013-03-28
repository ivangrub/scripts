#!/usr/bin/env python

import pysam as pys
import argparse
import os

def gen2bin(read,chrlist,headers,bin,chrt):
	s = read.rname
	chrom = chrt[s]
	ind = int(round(read.pos/float(bin-1)))
	chrst = chrlist[chrom]
	if chrst+ind >= len(headers):
		x = headers[len(headers)-1]['SN'].split('!')
		y = len(headers)-1
	else:
		x = headers[chrst+ind]['SN'].split('!')
		y = chrst+ind
	
	
	off = read.pos - int(x[2])+1
	k = 1
	while chrom != x[1] or off < 0:
		#print read.pos,x
		x = headers[y-k]['SN'].split('!')
		off = read.pos - int(x[2]) + 1
		k += 1

	if off > int(headers[y-k]['LN']) or off < 0:
		print 'ERROR',read

	read.pos = off
	read.tid = y-k+1
	return read	
	
	
def head2chr(header):
	chr = {}
	i = 0
	j = 0
	for line in header:
		x = line['SN']
		s = x.strip().split('!')
		if i == 0:
			chrom = s[1]
			chr[chrom] = i
			i += 1
		else:
			if s[1] == chrom:
				i += 1
				j += 1
			else:
				chrom = s[1]
				chr[chrom] = i
				i += 1
				j = 1
	chr[chrom] = i - j
	return chr

def offset(header):
	i = 0
	for line in header:	
		x = line['SN']
		s = x.strip().split('!')
		if i == 0:
			start = int(s[2])
			i = 1
		else:
			off = int(s[2]) - start + 1
			return off

parser = argparse.ArgumentParser(description='Put mapped reads into BAM format in chromosomal and express coordinates')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = '')
parser.add_argument('-l',default='50')
parser.add_argument('-g',default='mm9_norandom')
parser.add_argument('-o',default='outputted')
parser.add_argument('-b',default=1000)
args = parser.parse_args()

path = os.environ['EXPRESS_FILES']
if args.r is '-':
	infile = pys.Samfile('-','r')
else:
	infile = pys.Samfile(args.r,'rb')

conv = pys.Samfile('%s/Header_%s_%s_%s.sam' % (path,args.g,args.b,args.l),'r')
if args.r is '-':
	outbam = pys.Samfile('%s.bam' % args.o,'wb',template = infile)
convbam = pys.Samfile('%s_converted.bam' % args.o,'wb',template = conv)

headerlist = conv.header['SQ']
chr_tuple = infile.references
chrindex = head2chr(headerlist)
binning = offset(headerlist)

for line in infile:
	if line.is_unmapped:
		continue
	if args.r is '-':
		outbam.write(line)
	readmod = gen2bin(line,chrindex,headerlist,binning,chr_tuple)
	convbam.write(readmod)

if args.r is '-':
	outbam.close()
convbam.close()
infile.close()
conv.close()


