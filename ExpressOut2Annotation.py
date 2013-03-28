#!/usr/bin/env python

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Create a genome fasta file that will be used to create the bowtie index and used as an input eXpress.')
parser.add_argument('-g',help = 'Enter the name of the fasta file. The default is ref_mm9.fa',default = 'ref_mm9.fa')
args = parser.parse_args()


file = open('results.xprs')
fileout = open('results_mod.xprs','w')
fasta = SeqIO.parse('/Users/ivang/Bioinformatics/bowtie-index/%s' % args.g,'fasta')

fastaid = []
for line in fasta:
	s = line.description.strip().split()
	fastaid.append(s[2])

i = 0
for line in file:
	if i == 0:
		i += 1
		fileout.write(line)
		continue
	s = line.strip().split()
	s[1] = fastaid[int(s[1])-1]
	s.append('\n')
	x ='\t'.join(s)
	fileout.write(x)

fasta.close()
file.close()
fileout.close()