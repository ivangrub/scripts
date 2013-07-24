#!/usr/bin/env python

annots = ["refgene","ensgene","cufflinks_denovo"]


for type in annots:
	gene_annot = open("GeneChip_%s.txt" % type,'r')
	out = open("GeneChIP_%s_wRPKM.txt" % type,'w')
	if type == 'cufflinks_denovo':
		type = 'denovo'

	id = {}
	file = open("express_%s/results_mod.xprs" % type,'r')
	i = 0
	for line in file:
		if i == 0:
			i += 1
			continue
		s = line.strip().split()
		id[s[1]] = s[10]

	i = 0
	for line in gene_annot:
		if i == 0:
			i += 1
			s = line.strip().split()
			x = "\t".join(s)
			header = x+"\t"+"RPKM"
			out.write("%s\n" % header)
			continue
		s = line.strip().split()
		print s
		try:
			rpkm = id[s[0]]
		except KeyError:
			rpkm = "NaN"
		x = "\t".join(s) + "\t" + rpkm
		out.write("%s\n" % x)

	out.close()
	gene_annot.close()



	
