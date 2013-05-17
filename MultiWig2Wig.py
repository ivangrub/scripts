#!/usr/bin/env python

x = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY','chrM']

outfile = open('Combinedwiggle.wig','w')

j = 0
for i in x:
	k = 0
	wig = open('%s.wig' % i)
	for line in wig
		if  j== 0 and k == 0:
			outfile.write(line)
			j += 1
			k += 1
			continue
		elif k == 0:
			k += 1
			continue

		outfile.write(line)

	wig.close()

outfile.close()