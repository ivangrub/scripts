#!/usr/bin/env python

from Bio import SeqIO


out = open('RTCC07_WT9.trimmed.fastq','w')

for read in SeqIO.parse('RTCC07_WT9.fastq','fastq'):
	s = read.seq 
	edit.seq = s[0:49]
	SeqIO.write(edit.seq,out,'fastq')

out.close()


