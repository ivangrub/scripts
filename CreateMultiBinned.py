#!/usr/bin/env python

import os

bins = [500,1000,10000]

for i in bins:
	overlap = [36,i*.25,i*.5,i*.75]
	for j in overlap:
		os.system('./BinnedFasta.py -g Histone_repeat_insert3.fa -l %d -b %d' % (j,i))