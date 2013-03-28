#!/usr/bin/env python

import os

insert = ['Histone_repeat_insert1','Histone_repeat_insert2','Histone_repeat_insert3']
for k in insert:	
	os.system('bowtie-build -f %s.fa %s' % (k,k))