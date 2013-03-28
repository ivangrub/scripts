#!/usr/bin/env python

bin = 50

msize = open('mm9_norandom.size')

chr = {}
for line in msize:
	s = line.strip().split()
	chr[s[0]] = s[1]

bed = open('%dbin.bed' % bin,'w')

i = 0
for j in chr.keys():
	k = 0
	while k <= chr[j]:
		bed.write('%s\t%d\t%d\n' %(j,k,k+bin-1))
		k += bin -1

msize.close()
bed.close()
