#!/usr/bin/env python

from Bio import SeqIO
import pysam as pys
import numpy as np
import argparse

def build_fasta():
	"""Initialize a dictionary of each reference and the relevant sequence"""
	chr = {}
	for x in SeqIO.parse(args.g,'fasta'):
		chr[x.id] = x.seq
	return chr 

def MakeSAMHeader(chrom,gen):
	"""Print a SAM header that will be used in the final SAM file that is outputted after the projection of the reads back onto the genome"""
	ref = []
	header = {}
	header['HD'] = {'VN':'1.0'}
	for i in chrom.keys():
		ref.append({'LN':len(chrom[i]), 'SN':i})
	header['SQ'] = ref
	samheader = pys.Samfile('PostExpress_Headers_%s.sam' % gen,'wh',header = header)
	samheader.close()
	
def ReadLength(file):
	"""Determine the length of read outputted from the sequencing"""
	for x in SeqIO.parse(file,'fastq'):
		return len(x.seq)
	
def print_fasta(pos,seq):
	"""Print the new FASTA file that will be used to build the bowtie index and be the target sequence in eXpress"""
	binnum = pos.split('!')
	chr = binnum[1]
	coor1 = int(binnum[2])-1
	coor2 = int(binnum[3])-1
	NEWFasta.write('>%s\n%s\n\n' % (pos, seq[coor1:coor2].upper()))

parser = argparse.ArgumentParser(description='Create a genome fasta file that will be used to create the bowtie index and used as an input eXpress.')
parser.add_argument('-g',help = 'Enter the name of the fasta file. The default is mm9.fa',default = 'mm9.fa')
parser.add_argument('-f',help = 'Input the fastq file that will be used by bowtie.',default = None)
parser.add_argument('-l',help = 'Input the length of the fasta read',default = None)
args = parser.parse_args()

if args.f is None and args.l is None:
	error("You need to input a fastq file or at least the length of the read")
	
# Start of the main function
genome = args.g.strip().split('.')
NEWFasta = open('eXpress_%s.fa' % genome[0],'w')

# Read in information about the chromosomes (Both sequence and length parameters)
if args.f is None:
	length = int(args.l)
else:
	length = ReadLength(args.f)

# Output the reference dictionary into chr and print a new SAM header that will be used to make the final SAM file
chr = build_fasta()
MakeSAMHeader(chr,genome[0])

# Bin the sequences into smaller fasta format sequences which will be used to create the bowtie index
k = 0
for key in chr.keys():
	left = np.arange(1,len(chr[key]),1000-length)
	right = np.arange(1000,len(chr[key])+1000,1000-length)
	
	print 'On %s' % key
	for j in range(len(left)):
		bin = 'bin%d!%s!%d!%d' % (k,key,left[j],right[j])
		print_fasta(bin,chr[key])
		k += 1

NEWFasta.close()			
											
