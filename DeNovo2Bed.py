#!/usr/bin/env python

from Bio import SeqIO

fasta = SeqIO.parse('cufflinks.denovo.fa','fasta')
bed = open('cufflinks.denovo.bed','w')

for line in fasta:
	s = line.description.strip().split()
	csplit = s[1].split(':')
	chr = csplit[0]
	coord = csplit[1].split('-')
	bed.write('%s\t%s\t%s\t%s\t%d\t%s\n' % (chr,coord[0],coord[1],s[2],0,'+'))

fasta.close()
bed.close()