#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Project SAM/BAM format reads back onto the appropriate genome and call significant peaks.')
parser.add_argument('-m',help ='DECOD Motif Output',default = '')
args = parser.parse_args()

name = open(args.m,'r')
outfile = open('%s.stampout' % args.m,'w')
for line in name:
	s = line.strip()
	if s.find('>Motif') >= 0 or s.find('A') == 0 or s.find('T') == 0 or s.find('G') == 0 or s.find('C') == 0:
		outfile.write('%s\n' % s)

outfile.close()
name.close()
