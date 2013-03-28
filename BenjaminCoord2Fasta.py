#!/usr/bin/env python

from Bio import SeqIO

fa = open('/Users/ivang/Bioinformatics/bowtie-index/dm3.fa','r')
coord = open('/Users/ivang/Downloads/databasesdivergentpromoters/bmc.txt','rU')
coord2fa = open('/Users/ivang/Desktop/DivergentPromoters.fasta','w')



seq = SeqIO.to_dict(SeqIO.parse(fa,'fasta'))

i = 0
for line in coord:
	if i == 0:
		i += 1
		continue
	s = line.strip().split()
	if int(s[10]) <= 0:
		continue
	chr = 'chr'+s[1]
	sequence = seq[chr].seq
	coord2fa.write('>%s %s:%s-%s\n%s\n' %(s[2],chr,s[4],s[7],sequence[int(s[4])-200:int(s[7])+200]))

fa.close()
coord.close()
coord2fa.close()
