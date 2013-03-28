#!/usr/bin/env python

exo = open('TFIIH.peaks.gff.txt','r')
exout = open('TFIIH.peaks.bed','w')

i = 0
for line in exo:
	if i == 0:
		i += 1
		continue
	s = line.strip().split('\t')
	x = '\t'.join([s[0],s[6],s[3],s[4],s[5]])
	exout.write('%s\n' % x)


