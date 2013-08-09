#!/usr/bin/env python

results = open('Express_TBP.results.xprs')
out = open('Express_TBP.bins.truth.bed','w')

i = 0
for line in results:
	if i == 0:
		i += 1
		continue
	s = line.strip().split()
	bin = s[1].split('!')
	try:
		out.write('%s\t%s\t%s\t%s\n' % (bin[1],bin[2],bin[3],s[13]))
	except IndexError:
		print line

out.close()
results.close()