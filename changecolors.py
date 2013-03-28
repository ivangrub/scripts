#!/usr/bin/env python

x = raw_input('What is the name of the file? Include the extension\n')

op = open(x,'r')
out = open('newtest.txt','w')
i = 0
for line in op:
	m = divmod(i,3)
	i += 1
	if m[1] == 0:
		out.write('\033[31m%s\033[0m' % line) # print in red
		continue
	if m[1] == 1:
		print '\033[34m%s\033[0m' % line # print in blue
		continue
	if m[1] == 2:
		print '\033[32m%s\033[0m' % line # print in green
		continue

op.close()
out.close()