#!/usr/bin/env python

import numpy as np

IP = ['CSEM_Pol2','CSEM_PI','CSEM_TBP','Unique_Pol2','Unique_PI','Unique_TBP','Express_Pol2','Express_PI','Express_TBP','NoExpress_Pol2','NoExpress_PI','NoExpress_TBP']

for i in IP:
	file = open('%s.peaks.bed' % i,'r')
	newfile = open('%s.broad.bed' % i,'w')
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
				enr = np.amax(enrich)
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

