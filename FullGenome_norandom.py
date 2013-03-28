#!/usr/bin/env python

from Bio import SeqIO

fa = open('/Users/ivang/Bioinformatics/bowtie-index/hg19.fa','r')
fanorandom = open('/Users/ivang/Bioinformatics/bowtie-index/hg19_norandom.fa','w')

for seq in SeqIO.parse(fa,'fasta'):
	if '_' in seq.id:
		continue
	else:
		fanorandom.write('>%s\n%s\n'% (seq.id,seq.seq))

fa.close()
fanorandom.close()