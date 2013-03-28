#!/usr/bin/env python

IP = 'Rad23b'

for i in IP:
	newname = open('%s_SPP_NarrowPeaks.txt','r')
	out = open('%s_Peaks.txt','w')
	
	for line in newname:
		s = line.strip().split()
		chr = '%s:%d-%d' % (s[0],s[1],s[2])
		out.write(chr)
	
	newname.close()
	out.close()