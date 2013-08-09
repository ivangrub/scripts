#!/usr/bin/env python

init = ["Express","CSEM","Unique"]
types = ["Express_CSEM","Express_Unique","CSEM_Unique","Express_CSEM_Unique"]
outfile = open('OverlapCounts.TBP.txt','w')

for first in init:
	text = open("/Users/ivang/Desktop/ExpressFilter/Peaks/Filtered/%s_TBP.filtered.bed" % first)

	i = 0
	for line in text:
		i += 1

	outfile.write('%s\t%d\n' % (first,i))
	text.close()

for typ in types:
	text = open("%s.TBP.bed" % typ)

	i = 0
	for line in text:
		i += 1

	outfile.write('%s\t%d\n' % (typ,i))
	text.close()

outfile.close()
