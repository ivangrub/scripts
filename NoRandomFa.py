#!/usr/bin/env python

from Bio import SeqIO

x = open('mm9.fa','r')
out = open('mm9_norandom.fa','w')

for chr in SeqIO.parse(x,'fasta'):
	if 'random' not in chr.id:
		out.write('>%s\n%s\n' % (chr.id,chr.seq))
		
out.close()
x.close()