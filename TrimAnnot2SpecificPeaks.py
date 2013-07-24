#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Project SAM/BAM format reads back onto the appropriate genome and call significant peaks.')
parser.add_argument('-p',help ='Read the peak file.',default = '')
parser.add_argument('-a',help = 'The annotation file',default = '')
parser.add_argument('-o',help = 'The updated annotation file for a subset of the peaks',default = 'output.txt')
args = parser.parse_args()

peaks = open(args.p)
annot = open(args.a)
annot_new = open(args.o,'w')

list = set()
for line in peaks:
	s = line.strip().split()
	name = '%s:%s-%s' % (s[0],s[1],s[2])
	list.append(name)

i = 0
for line in annot:
	if i == 0:
		annot_new.write(line)
		i += 1
		continue
	s = line.strip().split()
	name = '%s:%s-%s' % (s[0],s[1],s[2])
	if name in list:
		annot_new.write(line)

peaks.close()
annot.close()
annot_new.close()
