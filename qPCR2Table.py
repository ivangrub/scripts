#!/usr/bin/env python

import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Reorganize qPCR output file')
parser.add_argument('-r',help ='File name',default = None)
args = parser.parse_args()

if args.r is None:
	error("Input a file name")


file = open(args.r,'r')
new = open('%s.txt' % args.r[:-4],'w')

new.write('Primer\tCondition\tCt1\tCt2\tCt3\tAverage Ct\tSt Dev Ct\n')

k = 0
i = 1
ct = []
for line in file:
	s = line.strip().split(',')
	if k == 0: 
		if s[0] == 'Well':
			k = 1
		continue
	if s[4] != 'Undetermined':
		ct.append(float(s[4]))
	else:
		ct.append(np.nan)
	if i % 3 == 0:
		mean = np.mean(ct)
		st = np.std(ct)
		new.write('%s\t%s\t%f\t%f\t%f\t%f\t%f\n' % (s[1],s[3],ct[0],ct[1],ct[2],mean,st))
		ct = []
	i += 1
	
new.close()
file.close()