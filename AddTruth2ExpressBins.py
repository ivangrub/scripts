#!/usr/bin/env python

results = open('results.xprs')
out = open('CoordwTruth.txt','w')

i = 0
for line in results:
	if i == 0:
		out.write('Chr\tCoord1\tCoord2\tSolvable\n')
		i += 1
		continue
	s = line.strip().split()
	bin = s[1].split('!')
	out.write('%s\t%s\t%s\t%s\n' % (bin[1],bin[2],bin[3],s[13])

out.close()
results.close()