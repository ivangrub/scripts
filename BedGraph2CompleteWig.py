#!/usr/bin/env python

bed = open('chip.bedgraph','r')
out = open('chip.complete.wig','w')
i = 0

for line in bed:
	if i == 0:
		i += 1
		continue

	if i == 1:
		s = line.strip().split()
		bin = int(s[1])
		y = [s[0],str(int(s[1])-bin),s[1],s[2]]
		out.write('\t'.join(y))
		i += 1
		continue

	s = line.strip().split()
	y = [s[0],str(int(s[1])-bin),s[1],s[2]]
	out.write('\t'.join(y))

out.close()
bed.close()