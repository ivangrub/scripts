#!/usr/bin/env python

import pysam as pys
import numpy as np
import argparse

def GenomeLen():
	"""Open the new SAM header file and create a binned array of zeros for each chromosome"""
	
	global sep
	genome = {}
	
	samfile = pys.Samfile('/Users/ivang/Bioinformatics/Bioinformatics/express/Post-eXpress_Header_%s.sam' % args.g,'r')
	header = samfile.header['SQ']
	
	# Iterate through the reference sequences in the header
	for i in header:
		genome[i['SN']] = np.zeros(np.floor(i['LN']/int(args.w))+sep)
	
	samfile.close()
	return genome,header

def AddtoBin(read):
	global chip, sep, chipF, chipR
	
	# Check the strand
	strand = read[1]
	
	# Figure out which reference sequence the read mapped to		
	bin = read[0].split('!')
	
	chr = bin[0]
	start = int(bin[1])
		
	st = int(np.floor(start/int(args.w)))
		
	# Create a coordinate index to update the initial zeros vector array in the chip dictionary
	
	if strand == '+':
		pp = np.arange(st,st+sep-1)
	else:
		pp = np.arange(st-sep+1+int(np.floor(36/int(args.w))),st+int(np.floor(36/int(args.w))))
		
	# Check boundary conditions
	if min(pp) > 1 and max(pp) < len(chip[chr]):
		chip[chr][pp] += 1
		if strand == '+':
			chipF[chr][pp] += 1
		else:
			chipR[chr][pp] += 1
		
def NormalizeReads(nreads):
	global chip, chipF, chipR
	# Normalize to 10 thousand reads
	ratio = nreads/1e4
	for i in chip.keys():
		chip[i] = chip[i]/ratio	
		chipF[i] = chipF[i]/ratio
		chipR[i] = chipR[i]/ratio

		
def ReadReads():
	"""Pass the reads that fit the predefined region"""	
	global region
	
	counts = open('%s.Counts.txt' % args.r,'r')
	locations = open('%s.ReadsAndLocations.txt' % args.r,'r')
			
	# Create a tuple of the reference sequences in the file
	print 'Making the location dictionary'
	DictofLocations = MakeReadDictionary(locations)
	nreads = 0
	print 'Starting to go through all of the reads in the region'	
	for line in counts:
		s = line.strip().split('\t')
		if s[0] == 'ID':
			continue
		if region == 'OnlyinRegion':
			if float(s[3])/float(s[1]) == 1 and int(s[2]) == 0 and int(s[1]) > 10 and int(s[1]) < 21:
				allloc = DictofLocations[s[0]]
				for i in range(len(allloc)):
					AddtoBin(allloc[i])
		else:
			if float(s[3])/float(s[1]) < 1 and int(s[2]) == 0 and int(s[1]) > 10 and int(s[1]) < 21:
				allloc = DictofLocations[s[0]]
				for i in range(len(allloc)):
					AddtoBin(allloc[i])
					
		nreads += 1
		if nreads % 1e4 == 0:
			print 'Processed %d reads' % nreads
			
	NormalizeReads(nreads)
	counts.close()
	locations.close()

def ParseHeader(hd):
	"""Parse the header list of dictionaries"""
	chr_list = []
	chr_len = []
	for i in hd:
		chr_list.append(i['SN'])
		chr_len.append(i['LN'])
	return chr_list,chr_len	
		
def PrintWiggle(reads,name):
	"""Print a wiggle file of the projected reads so that it can be visualized on the UCSC genome browser"""
	global chrlist,chrlen,region
	
	print 'Printing the %s wiggle file' % name
	f = open('%s.%s-%s.%s-11.wig' % (name,args.o,args.reg,args.g),'w')
	f.write('track type=bedGraph name="%s-%s_%s-11" description="%s-%s_%s-11" visibility=full graphType=bar\n' % (args.o,region,name,args.o,region,name))
	for i in reads.keys():
		length = chrlen[chrlist.index(i)]
		lastpos = np.nan
		lastval = np.nan
		X = np.arange(int(args.w),int(args.w)*(len(reads[i])-17),int(args.w))
		dat = reads[i]
		for zz in range(len(X)):
			d = abs(dat[zz] - lastval)
			if not np.isfinite(d) or d > 1e-1:
				if np.isfinite(lastval):
					if X[zz] > length:
						continue
					f.write('%s %d %d %.1f\n' % (i, lastpos, X[zz], lastval))
				lastval = dat[zz]
				lastpos = X[zz]
	f.close()

def MakeReadDictionary(loc):
	reads = {}
	for line in loc:
		s = line.strip().split('\t')
		if s[0] == 'ID':
			continue
		try:
			reads[s[0]].append(['%s!%s' % (s[1],s[2]),s[4]])
		except KeyError:
			reads[s[0]] = [['%s!%s' % (s[1],s[2]),s[4]]]
	
	return reads	
	 	
parser = argparse.ArgumentParser(description='Project SAM/BAM format reads back onto the appropriate genome and call significant peaks.')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = '')
parser.add_argument('-g',help = 'The genome abbreviation. Mouse is mm9. Human is hg18. Drosophila is dm3. The default is mm9.',default = 'mm9')
parser.add_argument('-w',help = 'The base pair binning. The default = 50bp.',default = 50)
parser.add_argument('-d',help = 'The DNA fragment length. The default = 250bp.',default = 250)
parser.add_argument('-o',help = 'The prefix of the output files.',default = 'test')
parser.add_argument('-reg',help = 'The prefix of the output files.',default = 'only')

args = parser.parse_args()
if args.reg == 'only':
	region = 'OnlyinRegion'
else:
	region = 'ElsewhereToo'
	
# Binning
sep = int(np.ceil(int(args.d)/int(args.w)))

# Allow for the summed input along with specific strands
v = ['chip','chipF','chipR']
chip,samheader = GenomeLen()
chipF, ign = GenomeLen()
chipR, ign = GenomeLen()
chrlist,chrlen = ParseHeader(samheader)

ReadReads()
for i in v:
	c = eval(i)
	PrintWiggle(c,i)
