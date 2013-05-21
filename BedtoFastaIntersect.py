#!/usr/bin/env python

from Bio import SeqIO

peaks = open('96h_T7_and_MEF_TAF7.peaks.bed','r')
seq = open('96h_T7_and_MEF_TAF7.peaks.txt','r')

newfasta = open('file1.filtered.fa','w')

record_dict = SeqIO.to_dict(SeqIO.parse(seq, "fasta"))

for line in peaks:
	s = line.strip().split('\t')
	id = '%s:%s-%s' % (s[0],s[1],s[2])
	newfasta.write('>%s\n%s\n' %(id,record_dict[id].seq))

peaks.close()
seq.close()
newfasta.close()