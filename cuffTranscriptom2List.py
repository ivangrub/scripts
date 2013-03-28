#!/usr/bin/env python

from Bio import SeqIO

cuff = open('cuffTranscriptome.fa','r')
out = open('CuffTranscriptIDs','w')
for record in SeqIO.parse(cuff,'fasta'):
	s = record.description.strip().split()
	out.write('%s\t%s\t%s\n' % (s[0],s[1],s[2]))

out.close()
cuff.close()
