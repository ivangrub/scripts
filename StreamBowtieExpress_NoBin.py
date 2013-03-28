#!/usr/bin/env python

import os

insert = ['Histone_repeat_insert1','Histone_repeat_insert2','Histone_repeat_insert3']

print 'Mapping with bowtie'
for k in insert:
	if os.path.exists('/Users/ivang/Desktop/NoExpress/GRO_%s.bam' % k):
		continue
	os.system('bzip2 -dc /Volumes/Ivan/ChIP-Seq/Benjamin/SRR073008.fastq.bz2 | bowtie -c -q -v 2 -aS %s - | samtools view -bS - > ~/Desktop/NoExpress/GRO_%s.bam' % (k,k))
			
#print 'Running eXpress'			
#for k in insert:
#	os.system('express -o ~/Desktop/NoExpress/GRO_%s --output-align-prob %s.fa ~/Desktop/NoExpress/GRO_%s.bam' % (k,k,k))

print 'Creating wiggles'			
for k in insert:
	if os.path.exists('/Users/ivang/Desktop/NoExpress/chip.NoeXpress.%s.wig' % k):
		continue
	os.system('./express2wiggle.py -r ~/Desktop/NoExpress/GRO_%s.bam -g %s -w 10 -d 36 -o NoeXpress_%s ' % (k,k,k))
	os.system('rm ~/Desktop/NoExpress/*pkl')
	
 