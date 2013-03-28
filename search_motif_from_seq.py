#!/usr/bin/env python
'''searhing motif represented by IUPAC'''

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
try:
	import motility as mt
except:
	print >>sys.stderr, "Error: Cannot find \'motility\' module. Download and install 'motility' from \'http://cartwheel.caltech.edu/motility/\'"

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from bx.binned_array import BinnedArray
from bx_extras.fpconst import isNaN
from bx.bitset_utils import *
from bx.bbi.bigwig_file import BigWigFile

#built in modules
from qcmodule import fasta
from qcmodule import cigar
from qcmodule import BED
from qcmodule import motif
#changes to the paths

__author__ = "Liguo Wang"
__copyright__ = "Copyright 2012, Wei Li's Lab"
__credits__ = []
__license__ = "GPL"
__version__ = "2.0"
__maintainer__ = "Liguo Wang"
__email__ = "liguow@bcm.edu"
__status__ = "Development" #Prototype or Production



def main():
	usage="%prog [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-s","--sequence",action="store",type="string",dest="dna_seq",help="DNA sequence in fasta format")
	parser.add_option("-m","--motif",action="store",type="string",dest="motif_iupac",help="Tab (space) separated two-column file: Motif_ID <space> Motif_IUPAC. Motif_ID should be unique. For example:Motif.7.4	TGTWCHH\nMotif.6.3	RGWACA\nMotif.6.2	TGTWCW")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output",help="Output file")
	parser.add_option("-n","--mismatch",action="store",type="int",dest="mismatch_num",default=0, help="Number of mismaatch. default=%default")

	(options,args)=parser.parse_args()

	if not (options.dna_seq and options.output and options.motif_iupac):
		parser.print_help()
		sys.exit(0)
	FOUT1 = open(options.output + '.seq2motif.xls','w')
	FOUT2 = open(options.output + '.motif2seq.xls','w')
	motifs={}
	
	print >>sys.stderr, "Reading motif file " + options.motif_iupac
	for line in open(options.motif_iupac,'r'):
		line=line.strip()
		if line.startswith(('#',' ','\n')):continue
		id,iupac = line.split()
		motifs[id] = iupac
	print >>sys.stderr, "Search motifs for each sequence "	
	for line in open(options.dna_seq,'r'):
		line=line.strip()
		if line.startswith(('#',' ','\n')):continue
		if line.startswith('>'):
			print >>FOUT1, line[1:] + '\t',
			continue
		print >>FOUT1, line + '\t',
		for motif in sorted(motifs):
			if len(motifs[motif]) > len(line):
				continue
			found = mt.find_iupac(line, motifs[motif], options.mismatch_num)
			print >>FOUT1, motif + ";" + motifs[motif] + ';' + str(len(found)) + '\t',
		print >>FOUT1
	
	print >>sys.stderr, "Search sequences for each motif"
	
	for motif in sorted(motifs):
		count = 0
		SEQ = open(options.dna_seq,'r')
		for line in SEQ:
			line=line.strip()
			if line.startswith(('#',' ','\n','>')):continue
			if len(motifs[motif]) > len(line): continue
			found = mt.find_iupac(line, motifs[motif], options.mismatch_num)
			if len(found)>0:count +=1
		SEQ.close()
		print >>FOUT2, motif + '\t' + motifs[motif]  + '\t' + str(count)
	
	
			
										
if __name__ == '__main__':
	main()
