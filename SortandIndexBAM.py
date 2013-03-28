#!/usr/bin/env python

import os

IP = ['CC-PI','CC-Rad23b','ES-PI','ES-TAF7','ES-TBP']

for i in IP:
	os.system('samtools sort Mostlikely_%s.mm9_norandom.ml.bam %s.ml.sorted' % (i,i))
	os.system('samtools index %s.ml.sorted.bam' % i)
	os.system('igvtools count -z 10 -w 50 -e 250 %s.ml.sorted.bam %s.ml.sorted.tdf mm9' % (i,i))
	