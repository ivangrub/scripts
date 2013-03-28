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
	parser.add_option("-b","--bed",action="store",type="string",dest="bed_file",help="bed format file")
	parser.add_option("-r","--reference",action="store",type="string",dest="fasta_file",help="reference sequece in fasta format")
	parser.add_option("-m","--motif",action="store",type="string",dest="motif_iupac",help="Tab (space) separated two-column file: Motif_ID <space> Motif_IUPAC. Motif_ID should be unique. For example:Motif.7.4	TGTWCHH\nMotif.6.3	RGWACA\nMotif.6.2	TGTWCW")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output",help="Output file")
	parser.add_option("-n","--mismatch",action="store",type="int",dest="mismatch_num",default=0, help="Number of mismaatch. default=%default")

	(options,args)=parser.parse_args()

	if not (options.bed_file and options.output and options.motif_iupac and options.fasta_file):
		parser.print_help()
		sys.exit(0)
	
	obj = fasta.Fasta(options.fasta_file)
	FOUT1 = open(options.output,'w')
	
	for line in open(options.bed_file,'r'):
		line=line.strip()
		fields=line.split()
		chrom = fields[0]
		st = int(fields[1])
		end = int(fields[2])	
		motif = obj.search_motif_iupac(chrom, st, end, options.motif_iupac,options.mismatch_num,1)
		if len(motif) >0:
			print >>FOUT1, line + '\t' + motif[0] 	
		else:
			print >>FOUT1, line + '\tNA'
			
										
if __name__ == '__main__':
	main()
