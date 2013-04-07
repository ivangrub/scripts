#!/usr/bin/env python

import pysam as pys
import argparse
import os

def gen2bin(infile,read,chrlist,headers,bin,headlist):
	chrom = infile.getrname(read.tid)
	ind = int(round(read.pos/float(bin-1)))
	chrst = chrlist[chrom]

	try:
		name = convbam.getrname(chrst+ind)
		x = name.split('!')
		y = chrst+ind
	except TypeError:
		length = convbam.nreferences
		name = convbam.getrname(length -1)
		x = name.split('!')
		y = length-1
	
	off = read.pos - int(x[2])+1
	k = 1
	
	while chrom != x[1] or off < 0:
		name = convbam.getrname(y-k)
		x = name.split('!')
		off = read.pos - int(x[2]) + 1
		k += 1

	if off > int(int(args.b) - 1) or off < 0:
		print 'ERROR',read

	read.pos = off
	read.tid = convbam.gettid(name)
	indice = [-1,0,1]
	for ran in indice:
		try:
			if (convbam.getrname(y-k+1+ran) not in headlist):
				headlist.add(convbam.getrname(y-k+1+ran))
			else:
				continue
		except TypeError:
			continue

	return read,headlist	

def MakeHeader(binheader,binsize):
	ref = []
	header = {}
	nh = {}
	header['HD'] = {'VN':'1.0'}
	index = 0
	for i in binheader:
		ref.append({'LN':int(binsize)-1,'SN':i})
		nh[i] = index
		index += 1
	header['SQ'] = ref
	return header,nh	

def TrimBAM(conv,oheader,nheader,translate):
	trimmed = pys.Samfile('%s_trimmed.bam' % args.o,'wb',header=nheader)
	for read in conv:
		binid = oheader[read.tid]["SN"]
		index = translate[binid]
		read.tid = index
		trimmed.write(read)
	trimmed.close()

def head2chr(header):
	print 'Getting header index'
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
convbam = pys.Samfile('%s_converted.bam' % args.o,'wb',header = conv)

bindict = {}
index = 0
for ref in conv.header['SQ']:
	bindict[ref["SN"]] = index
	index += 1

print 'Building conversion headers'
chrindex = head2chr(headerlist)

print 'Getting the offset'
binning = offset(headerlist)

print 'Binning the alignments'
binheader = set()
for line in infile:
	if line.is_unmapped:
		continue
	if args.r is '-':
		outbam.write(line)
	readmod,binheader = gen2bin(infile,line,chrindex,bindict,binning,binheader)
	convbam.write(readmod)

convbam.close()
conv.close()
conv = pys.Samfile('%s_converted.bam' % args.o,'rb')

print 'Making new header'
newheader, headdict = MakeHeader(binheader,args.b)

print 'Trimming and updating tid'
TrimBAM(conv,headerlist,newheader,headdict)
conv.close()

if args.r is '-':
	outbam.close()

infile.close()



