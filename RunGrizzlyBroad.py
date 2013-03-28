#!/usr/bin/env python

import os

dir = os.listdir('/Users/ivang/Desktop')

for j in dir:
	if 'chip.' in j and 'peaks.bed' in j and 'overlap' not in j:
		os.system('Grizzly_OnlyBroad.py -i %s' % j)