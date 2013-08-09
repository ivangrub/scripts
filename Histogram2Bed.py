#!/usr/bin/env python

old_hist = open('Express_TBP.histogram.txt')
new_hist = open('Express_TBP.histogram.bed','w')

for line in old_hist:
	s = line.strip().split()
	bed = s[0].split(':')
	new = "\t".join(bed)
	new_hist.write(new+"\t"+"\t".join(s[1:])+"\n")

old_hist.close()
new_hist.close()