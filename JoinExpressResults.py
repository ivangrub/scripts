#!/usr/bin/env python

file1 = '/Users/ivang/Desktop/HY_3.results.xprs'
file2 = '/Users/ivang/Desktop/HY_5.results.xprs'
s1 = open(file1)
s2 = open(file2)

fpkm = {}
i = 0
for line in s1:
	if i == 0:
		i += 1
		continue
	s = line.strip().split()
	fpkm[s[1]] = float(s[10])

i = 0
for line in s2:
	if i == 0:
		i += 1
		continue
	s = line.strip().split()
	fpkm[s[1]].append(float(s[10]))

s1.close()
s2.close()

out = open('HY_3and5.txt','w')

out.write('%s\t%s\t%s\n' % ('Genes',file1,file2))
for key in fpkm.keys():
	vals = fpkm[key]
	out.write('%s\t%f\t%f\n' % (key,vals[0],vals[1]))

out.close()