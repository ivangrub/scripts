#!/usr/bin/env python

# Convert Tommy Kaplan's Grizzly Peak (ChIP-Seq Peak calling program) from MATLAB into 
# python for integration with eXpress.

import numpy as np
import pysam as pys
import argparse
import glob
import pickle
import os

def InitBam():
	"""Open the new SAM header file and create a binned array of zeros for each chromosome"""
	
	genome = {}
	path = os.environ['EXPRESS_FILES']
	samfile = pys.Samfile('%s/Post-eXpress_Header_%s.sam' % (path,args.g),'r')
	header = samfile.header['SQ']
	
	# Iterate through the reference sequences in the header
	j = 0
	for i in header:
		genome[i["SN"]] = j
		j += 1
	
	samfile.close()
	return genome
		
def ReadReads(dir,chrlist):
	"""Pass in the BAM/SAM file reads and project them back onto the genome"""	
	path = os.environ['EXPRESS_FILES']
	samfile = pys.Samfile('%s/Post-eXpress_Header_%s.sam' % (path,args.g),'r')
	outfile = pys.Samfile('%s/Post-eXpress-%s.%s.bam' %(dir,args.o,args.g),'wb',template = samfile)
	
	if '.bam' in args.r:
		r = pys.Samfile(args.r,'rb')
	else:
		r = pys.Samfile(args.r,'r')
			
	refdict = r.references	
	for seq in r:	
		if seq.is_unmapped:
			continue
		bin = refdict[seq.tid]
		st = int(seq.pos)

		chrom = bin.split('!')
		seq.tid = chrlist[chrom[1]]
		seq.pos = st+int(chrom[2])-1
		outfile.write(seq)

		
	r.close()
	samfile.close()
	outfile.close() 
	

parser = argparse.ArgumentParser(description='Project SAM/BAM format reads back onto the appropriate genome and call significant peaks.')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = '')
parser.add_argument('-o',help = 'The prefix of the output files.',default = 'test')
args = parser.parse_args()

# Establish express output directory and print new files there
dir = args.r.split('/')
dirname = '/'.join(dir[:-1])

# Allow for the summed input along with specific strands

chrindex = InitBAM()

ReadReads(dirname,chrindex)

