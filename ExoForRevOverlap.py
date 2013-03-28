#!/usr/bin/env python

import numpy as np

forward = open('forward.bedgraph','r')
reverse = open('reverse.bedgraph','r')

fchr = []
fr = []
for line in forward:
	s = line.strip().split()
	fchr.append(s[0])
	fr.append(s[2])

rchr = []
rl = []
rr = []
for line in reverse:
	s = line.strip().split()
	rchr.append(s[0])
	rl.append(s[1])
	rr.append(s[2])

forward.close()
reverse.close()

