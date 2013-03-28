#!/usr/bin/env python

import os

bins = [500,1000,10000]
insert = ['Histone_repeat_insert1','Histone_repeat_insert2','Histone_repeat_insert3']
for k in insert:
	for i in bins:
		overlap = [36,i*.25,i*.5,i*.75]
		for j in overlap:
			name = 'eXpress_%dbp_%d.%s' % (i,j,k)
			os.system('bowtie-build -f %s.fa %s' % (name,name))