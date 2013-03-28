#!/usr/bin/env python

from Bio import SeqIO
import os

path = os.environ['BOWTIE_INDEXES']
x = '/mm9_norandom.fa'
out = open(path+'/mm9_norandom_size.txt','w')
for seq in SeqIO.parse(path+x,'fasta'):
	l = len(seq.seq)
	out.write('%s\t%d\n' % (seq.description,l))

out.close()