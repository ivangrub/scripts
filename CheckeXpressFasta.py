#!/usr/bin/env python

from Bio import SeqIO

f = '/Users/ivang/Bioinformatics/express/eXpress_500bp_50.mm9_norandom.fa'

for seq in SeqIO.parse(f,'fasta'):
	if len(seq.seq) != 499:
		print 'ERROR',len(seq.seq),seq.description

