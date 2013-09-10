#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Create a genome fasta file that will be used to create the bowtie index and used as an input eXpress.')
parser.add_argument('-p',help = 'Prefix of the wiggle files',default = None)
parser.add_argument('-o',help = 'Prefix of combined wiggle file',default = None)
args = parser.parse_args()

if args.p == None:
	print 'Provide a prefix file'

x = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY','chrM']

outfile = open('%s.wig' % args.o,'w')

j = 0
for i in x:
	k = 0
	prefix = args.p+'_'+i
	wig = open('%s.wig' % i)
	for line in wig:
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