#!/usr/bin/env python
'''searhing motif represented by IUPAC (input is bed)'''

#import built-in modules
import os,sys
import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
import sets
import numpy
import subprocess

__author__ = "Liguo Wang"
__copyright__ = "Copyright 2012, Wei Li's Lab"
__credits__ = []
__license__ = "GPL"
__version__ = "2.0"
__maintainer__ = "Liguo Wang"
__email__ = "liguow@bcm.edu"
__status__ = "Development" #Prototype or Production


seqcode = {"A":1,"C":2,"G":3,"T":4,"N":5}
code = []
seq_num = 0
def main():
	usage="%prog [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--dna",action="store",type="string",dest="dna_seq",help="DNA sequence in fasta format. sequnce length must be fixed")
	parser.add_option("-o","--out",action="store",type="string",dest="output",help="output prefix")
	(options,args)=parser.parse_args()
	if not (options.dna_seq and options.output ):
		parser.print_help()
		sys.exit(0)
	
	ROUT=open(options.output + '.r','w')
	for line in open(options.dna_seq):
		line = line.strip()
		if line.startswith('>'):continue
		line=line.upper()
		seqlen = len(line)
		for base in line:
			try:
				code.append(seqcode[base])
			except:
				continue
	print >>ROUT, "pdf('%s')" % (options.output + '.heatmap.pdf')
	print >>ROUT, "v=c(" + ','.join([str(i) for i in code]) + ')'
	print >>ROUT, "mat=matrix(v,byrow=T,ncol=%d)" % (seqlen)
	print >>ROUT, "heatmap(mat,Rowv=NA,Colv=NA,scale=c(\"none\"),col=c('green3','blue','orange','red'),labRow='')"
	print >>ROUT, "dev.off()"
	ROUT.close()

	subprocess.call("Rscript " + options.output + '.r',shell=True)

if __name__=='__main__':
	main()
