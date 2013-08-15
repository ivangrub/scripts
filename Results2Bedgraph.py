#!/usr/bin/env python

results = open('results.xprs')
out = open('/Users/grubisic/Desktop/Express_PI.results.bedgraph','w')

i = 0
for line in results:
	if i == 0:
		i += 1
		continue

	s = line.strip().split()
	chrom = s[1].split('!')
	out.write('%s\t%s\t%s\n' % (chrom[1],chrom[2],s[6]))

out.close()
results.close()