#!/usr/bin/env python

express = open('Pol2.express_40.dm3.unmapped.sam','r')
dm3 = open('Pol2.dm3.unmapped.sam','r')

exid = []
for line in express:
	s = line.strip().split('\t')
	y = s[0].split(' ')
	exid.append(y[0])

express.close()

dm3id = []	
for line in dm3:
	s = line.strip().split('\t')
	y = s[0].split(' ')
	dm3id.append(y[0])
	
dm3.close()

diffreads = open('DiffUnmapped.txt','w')

for x in dm3id:
	if x not in exid:
		diffreads.write('%s\n' % x)

diffreads.close()