#!/usr/bin/env python

import numpy as np

IP = ['AR','IgG','PolII','TBP','TAF7','TAF7L']

for i in IP:
	file = open('%s.peaks.bed' % i,'r')
	newfile = open('broad.%s.peaks.bed' % i,'w')
	j = 0
	enrich = []
	for line in file:
		if j == 0:
			j += 1
			continue
		s = line.strip().split(' ')
		
		try:
			x = int(s[3])
			if len(enrich) >= 1:
				enr = np.mean(enrich)
				y = [chrom,left,right,str(enr)]
				newfile.write('\t'.join(y)+'\n')
				enrich = []
			chrom = s[0]
			left = s[1]
			right = s[2]
		except ValueError:
			enrich.append(float(s[3]))

	file.close()
	newfile.close()

