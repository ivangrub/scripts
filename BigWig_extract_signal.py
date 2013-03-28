#!/usr/bin/env python

#import built-in modules
import os,sys
import re
import string
from optparse import OptionParser
import warnings
import string
import collections
from numpy import mean,median,std

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from bx.bbi.bigwig_file import BigWigFile

#import my own modules
from qcmodule import SAM
from qcmodule import fickett
from qcmodule import BED
from qcmodule import wiggle
from qcmodule import fasta
from qcmodule import orf
from qcmodule import quantile
#changing history to this module


def profile_bwfile(inbed,bwfile):
	'''retrieve signal from bigwig file for each entry in input bed file'''
	bw = BigWigFile( file=open( bwfile ) )
	
	for line in open(inbed):
		bw_signal=[]
		try:
			if line.startswith('#'):continue
			if line.startswith('track'):continue
			if line.startswith('browser'):continue
			if not line.strip():
				continue
			else:
				line = line.rstrip('\r\n')
				fields = line.split()
				chrom = fields[0]
				start = int(fields[1])
				end = int(fields[2])
		except:
			print >>sys.stderr,"Must be  chrom [space] start [space] end: " + line,
			continue		
		bw_signal.extend(bw.get_as_array(chrom,start,end))
		print chrom +'\t'+ str(start) +'\t'+ str(end) + '\t' + ','.join(str(i) for i in bw_signal)		


def main():
	if len(sys.argv) !=3:
		print >>sys.stderr, "\n================================================================================"
		print >>sys.stderr, "Usage: " + sys.argv[0] + "    input.bed    input.bw"
		print >>sys.stderr, " * Bed file: chrom  start  end"
		print >>sys.stderr, " * Direct output to STDOUT"
		print >>sys.stderr, "================================================================================"
		sys.exit(1)
	profile_bwfile(sys.argv[1],sys.argv[2])
if __name__ == '__main__':
	main()
