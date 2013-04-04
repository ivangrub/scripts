#!/usr/bin/env python

import pysam as pys
import argparse
import os

def gen2bin(read,chrlist,headers,bin,chrt,headlist):
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
		x = headers[y-k]['SN'].split('!')
		off = read.pos - int(x[2]) + 1
		k += 1

	if off > int(headers[y-k+1]['LN']) or off < 0:
		print 'ERROR',read

	read.pos = off
	read.tid = y-k+1
	indice = [-1,0,1]
	for ran in indice:
		try:
			headlist.append(headers[y-k+1+ran]['SN'])
		except IndexError:
			continue

	return read,headlist	

def MakeHeader(binheader,binsize):
	ref = []
	header = {}
	header['HD'] = {'VN':'1.0'}
	for i in binheader:
		ref.append({'LN':binsize, 'SN':i})
	header['SQ'] = ref
	return header	

def TrimBAM(conv,oheader,nheader):
	trimmed = pys.Samfile('%s_trimmed.bam' % args.o,'wb',template=nheader)
	for read in conv:
		binid = oheader[read.tid]
		index = nheader.index(binid)
		print read.tid,index
		read.tid = index
		trimmed.write(read)
	trimmed.close()

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

binheader = []
for line in infile:
	if line.is_unmapped:
		continue
	if args.r is '-':
		outbam.write(line)
	readmod,binheader = gen2bin(line,chrindex,headerlist,binning,chr_tuple,binheader)
	convbam.write(readmod)

convbam.close()
conv.close()
conv = pys.Samfile('%s_converted.bam' % args.o,'rb')
newheader = MakeHeader(binheader,binning)
print newheader
TrimBAM(conv,headerlist,newheader)
conv.close()

if args.r is '-':
	outbam.close()

infile.close()



