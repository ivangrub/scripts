#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Split out only broad peaks from Grizzly Peak peak file')
parser.add_argument('-i', type=str,help='Input Peak file in bed format')

args = parser.parse_args()

file = open(args.i,'r')
outfile = open('broad_%s' % args.i,'w')
narrow = open('narrow_%s' % args.i,'w')

j = 1
i = 0
for line in file:
	if i == 0:
		i += 1
		outfile.write(line)
		narrow.write(line)
		continue
	s = line.strip().split(' ')
	try:
		x = int(s[3])
		if x == j:
			outfile.write('\t'.join(s)+'\n')
			j += 1
	except ValueError:
		narrow.write('\t'.join(s)+'\n')


file.close()
outfile.close()
narrow.close()
		
	

