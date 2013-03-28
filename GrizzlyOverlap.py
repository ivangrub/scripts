#!/usr/bin/env python

import os

dir = os.listdir('/Users/ivang/Desktop')

peaks = []
for file in dir:
	if 'broad_chip.' in file and '.peaks.bed' in file:
		peaks.append(file)

for peak in peaks:
	for k in xrange(len(peaks)):
		if peak != peaks[k]:
			print peak, peaks[k]
			os.system('intersectBed -a ~/Desktop/%s -b ~/Desktop/%s > ~/Desktop/overlapped_%s_%s.bed' % (peak,peaks[k],peak[:-4],peaks[k][:-4]))