#!/usr/bin/env python

import os

overlap = {}
IP = ['TBP','TAF7','TAF7L','AR','PolII']
overlap['TAF7L'] = ['TBP','TAF7','AR']
overlap['TAF7'] = ['TBP','TAF7L','AR']
elem = ['promoter','proximal','distal']

dir = os.environ['BOWTIE_INDEXES']
for chip in IP:
	for cis in elem:
		os.system('intersectBed -a ~/Downloads/Grizzly_peaks/%s.filtered.peaks.bed -b %s/knownGene_%s.mm9.bed > ~/Downloads/Grizzly_peaks/%s.peaks_only_on.%s.bed' % (chip,dir,cis,chip,cis) )
		os.system('intersectBed -a %s/knownGene_%s.mm9.bed -b ~/Downloads/Grizzly_peaks/%s.filtered.peaks.bed > ~/Downloads/Grizzly_peaks/knownGene_Genes.%s_peaks_on.%s.bed' % (dir,cis,chip,chip,cis))
	os.system('intersectBed -a ~/Downloads/Grizzly_peaks/%s.peaks_on.proximal.bed -b ~/Downloads/Grizzly_peaks/%s.peaks_on.promoter.bed -v > ~/Downloads/Grizzly_peaks/%s.peaks_only_on.proximal.bed' %(chip,chip,chip))
	os.system('intersectBed -a ~/Downloads/Grizzly_peaks/%s.peaks_on.distal.bed -b ~/Downloads/Grizzly_peaks/%s.peaks_on.promoter.bed -v > ~/Downloads/Grizzly_peaks/%s.peaks_not_on_promoter.distal.bed' %(chip,chip,chip))
	os.system('intersectBed -a ~/Downloads/Grizzly_peaks/%s.peaks_not_on_promoter.distal.bed -b ~/Downloads/Grizzly_peaks/%s.peaks_on.proximal.bed -v > ~/Downloads/Grizzly_peaks/%s.peaks_only_on.distal.bed' %(chip,chip,chip))

for key in overlap.keys():
	for i in overlap[key]:
		os.system('intersectBed -a ~/Downloads/Grizzly_peaks/%s.filtered.peaks.bed -b ~/Downloads/Grizzly_peaks/%s.filtered.peaks.bed -wb > ~/Downloads/Grizzly_peaks/%s_on_%s.filtered.peaks.bed' % (key,i,key,i))
		os.system('intersectBed -a ~/Downloads/Grizzly_peaks/%s.filtered.peaks.bed -b %s/knownGene.mm9.bed -wb > ~/Downloads/Grizzly_peaks/%s_on_knownGene.filtered.peaks.bed' % (key,dir,key))
		for cis in elem:
			os.system('intersectBed -a ~/Downloads/Grizzly_peaks/%s.peaks_only_on.%s.bed -b ~/Downloads/Grizzly_peaks/%s.peaks_only_on.%s.bed -wb > ~/Downloads/Grizzly_peaks/%s_on_%s.on_%s.peaks.bed' % (key,cis,i,cis,key,i,cis))
