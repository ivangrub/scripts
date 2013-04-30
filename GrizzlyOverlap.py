#!/usr/bin/env python

import os

dir = os.listdir('/Volumes/Genomic1/IG_express/ES_Grizzly_Peaks')

peaks = []
for file in dir:
	if 'broad_chip.' in file and 'filtered.peaks.bed' in file:
		peaks.append(file)

for peak in peaks:
	for k in xrange(len(peaks)):
		if peak != peaks[k]:
			print peak, peaks[k]
			os.system('intersectBed -a /Volumes/Genomic1/IG_express/ES_Grizzly_Peaks/%s -b /Volumes/Genomic1/IG_express/ES_Grizzly_Peaks/%s > /Volumes/Genomic1/IG_express/ES_Grizzly_Peaks/overlapped_%s_%s.bed' % (peak,peaks[k],peak[:-4],peaks[k][:-4]))