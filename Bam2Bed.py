#!/usr/bin/env python

import pysam as pys
import argparse

parser = argparse.ArgumentParser(description='Project SAM/BAM format reads back onto the appropriate genome and call significant peaks.')
parser.add_argument('-b',help ='Read file name. BAM format.',default = '')

args = parser.parse_args()

x = pys.Samfile(args.b,'rb')
ref = x.references

j = 0
for seq in x:
	if j == 0:
		l = len(seq.seq)
		j += 1
	bin = ref[seq.rname]
	s = bin.split('!')
	chr = s[1]
	start = str(int(seq.pos) + int(s[2]))
	end = str(int(seq.pos) + int(s[2])+l-1)
	try:
		frac = str(seq.opt('XP'))
	except KeyError:
		frac = '1'
	print '\t'.join([chr,start,end,frac])
	
x.close()
	
	