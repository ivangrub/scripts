#!/usr/bin/env python

import os

bins = [1000,10000]
insert = ['Histone_repeat_insert1','Histone_repeat_insert2','Histone_repeat_insert3']

print 'Mapping with bowtie'
for k in insert:
	for i in bins:
		overlap = [36,i*.25,i*.5,i*.75]
		for j in overlap:
			name = 'eXpress_%dbp_%d.%s' % (i,j,k)
			if os.path.exists('/Users/ivang/Desktop/NewExpress/GRO_%s.bam' % name):
				continue
			os.system('bzip2 -dc /Volumes/Ivan/ChIP-Seq/Benjamin/SRR073008.fastq.bz2 | bowtie -c -q -v 2 -aS %s - | samtools view -bS - > ~/Desktop/NewExpress/GRO_%s.bam' % (name,name))
			
print 'Running eXpress'			
for k in insert:
	for i in bins:
		overlap = [36,i*.25,i*.5,i*.75]
		for j in overlap:
			name = 'eXpress_%dbp_%d.%s' % (i,j,k)
			if os.path.exists('/Users/ivang/Desktop/NewExpress/GRO_%s/hits.1.prob.bam' % name):
				continue
			os.system('express -o ~/Desktop/NewExpress/GRO_%s -B 1 --output-align-prob %s.fa ~/Desktop/NewExpress/GRO_%s.bam' % (name,name,name))

print 'Creating wiggles'			
for k in insert:
	for i in bins:
		overlap = [36,i*.25,i*.5,i*.75]
		for j in overlap:
			name = 'eXpress_%dbp_%d.%s' % (i,j,k)
			if os.path.exists('/Users/ivang/Desktop/NewExpress/GRO_%s/chip.PosteXpress_%s.%s.wig' % (name,name,k)):
				continue
			os.system('./express2wiggle.py -r ~/Desktop/NewExpress/GRO_%s/hits.1.prob.bam -g %s -w 10 -d 36 -o PosteXpress_%s ' % (name,k,name))
			os.system('rm ~/Desktop/NewExpress/GRO_%s/hits.1.prob.bam' % name)
 