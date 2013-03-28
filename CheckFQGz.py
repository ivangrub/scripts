#!/usr/bin/env python

import os

samples = ['HY_3','HY_4','HY_5']

for i in samples:
	print 'On sample %s' % i
	s = os.listdir('/Users/ivang/Downloads/Sample_%s_1_12' % i)
	for j in s:
		if '.fastq.gz' in j:
			print j
			os.system('gunzip /Users/ivang/Downloads/Sample_%s_1_12/%s' % (i,j))
	s = os.listdir('/Users/ivang/Downloads/Sample_%s_1_12' % i)
	for j in s:
		if '.gz' not in j and '.fastq' in j:
			os.system('cat /Users/ivang/Downloads/Sample_%s_1_12/%s > %s.fq' % (i,j,i))

	os.system('bzip2 -z /Users/ivang/Downloads/Sample_%s_1_12/%s.fq' % (i,i))

#for i in samples:
#	os.system('bzip2 -dc /Users/ivang/Downloads/Sample_%s/%s.fq | bowtie -c -q -v 2 -m 1 -S mm9_norandom - | samtools view -bS - > %s.mm9_norandom.bam' %(i,i,i))