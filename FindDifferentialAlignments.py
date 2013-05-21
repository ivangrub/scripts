#!/usr/bin/env python

pre = open('pre')
post = open('post')

diff = open('diff_alignments.txt','w')

predict = {}
for line in pre:
	s = line.strip().split('\t')
	predict[s[0]] = int(s[1])

for line in post:
	s = line.strip().split('\t')
	try:
		oldcount = predict[s[0]]
		if int(s[1]) != oldcount:
			diff.write('%s\t%s\t%s\n' % (s[0],str(oldcount),s[1]))
	except KeyError:
		diff.write('%s\t%d\t%s\n' % (s[0],0,s[1]))

pre.close()
post.close()
diff.close()

