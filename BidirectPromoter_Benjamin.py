#!/usr/bin/env python

import numpy as np


exp = open('supp.txt','rU')
prom = open('bmc.txt','rU')
out = open('Benjamin_Map.txt','w')

exp_map = {}
exp_list = []
i = 0
j = 0
for line in exp:
	s = line.strip().split('\t')
	try:
		if 'CG' not in s[4]:
			continue
	except IndexError:
		continue
	
	if i == 0:
		a= s[4]
		buf = a.split(':')
		id_gene = buf[0]
		exp_list.append(float(s[len(s)-1]))
		i += 1
		continue

	a = s[4]
	buf = a.split(':')
	if id_gene == buf[0]:
		exp_list.append(float(s[len(s)-1]))
	else:
		exp_map[id_gene] = [np.median(exp_list),max(exp_list),min(exp_list)]
		id_gene = buf[0]
		exp_list = []
		exp_list = [float(s[len(s)-1])]

for line in prom:
	s = line.strip().split('\t')
	if 'CG' not in s[2]:
		continue
	try:
		a = [s[2],str(exp_map[s[2]][0]),str(exp_map[s[2]][1]),str(exp_map[s[2]][2]),s[6],str(exp_map[s[6]][0]),str(exp_map[s[6]][1]),str(exp_map[s[6]][2]),str(s[len(s)-2]),'\n']
	except KeyError:
		continue
	x = '\t'.join(a)
	out.write(x)



exp.close()
prom.close()
out.close()