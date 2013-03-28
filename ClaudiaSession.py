#!/usr/bin/env python

from Bio import SeqIO
import os
import argparse

def load_fasta(refgen):
	# Calls the environmental variable BOWTIE_INDEXES
	path = os.environ['BOWTIE_INDEXES']

	# Create variable that defines the location and file name of reference genome
	fa = path + '/' + refgen
	#print refgen,path,fa

	# Initialize a dictionary to save the genome to memory
	genome = {}

	# Iterate over all sequences and save to genome
	for chr in SeqIO.parse(fa,'fasta'):
		genome[chr.id] = chr.seq

	#print genome.keys(),genome['chrM']
	return genome


# Get list of arguments
parser = argparse.ArgumentParser(description='Get sequences associated with peaks')
parser.add_argument('-p',help = 'Name of the peaks file',default = None)
parser.add_argument('-b',help = 'Total length of the sequence',default = 200)
parser.add_argument('-o',help = 'Name of the output fasta file',default = None)
parser.add_argument('-f',help = 'The reference fasta file',default = 'mm9_norandom.fa')
args = parser.parse_args()

if args.p is None or args.o is None:
	print "Tell me what file I need idiot"
	raise SystemExit

# Load sequence information from annotation into a dictionary called refgenome
refgenome = load_fasta(args.f)

# Open peaks file
peaks = open(args.p)
j = args.p.split('/')
dir = '/'.join(j[:-1])
peakseq = open(dir+'/'+args.o+'.fa','w')

bin = int(args.b)/2

i = 0
for line in peaks:
	if i == 0:
		i += 1 # same as i = i + 1
		continue
	pieces = line.strip().split('\t')
	chrom = pieces[0]
	left = int(pieces[1])-1
	right = int(pieces[2])-1
	mid = (right + left)/2
	sequence = refgenome[chrom][mid-100:mid+100]
	descript = '%s:%d-%d' % (chrom,left+1,right+1)
	peakseq.write('>%s\n%s\n' % (descript,sequence))


peaks.close()
peakseq.close()














