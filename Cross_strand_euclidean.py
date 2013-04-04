#!/usr/bin/env python
'''Calculate normalized Cross strand euclidean distance based on two bigWig files'''

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
import random
import numpy

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
from qcmodule import twoList
#changes to the paths


__author__ = "Liguo Wang"
__copyright__ = "Copyright 2012, Wei Li's Lab"
__credits__ = []
__license__ = "GPL"
__version__ = "2.0"
__maintainer__ = "Liguo Wang"
__email__ = "liguow@bcm.edu"
__status__ = "Development" #Prototype or Production

logo='''
   _____ _     ______          _____  
  / ____| |   |  ____|   /\   |  __ \ 
 | |    | |__ | |__     /  \  | |__) |
 | |    | '_ \|  __|   / /\ \ |  ___/ 
 | |____| | | | |____ / ____ \| |     
  \_____|_| |_|______/_/    \_\_|     
ChIP-Exo  Analysis  Pipeline (Cross Strand Distance)                                      
'''

def load_chromsize(file):
	'''read chrom.size file'''
	chromSize={}
	for line in open(file,'r'):
		if line.startswith('#'):continue
		if not line.strip():continue
		fields = line.strip().split()
		chromSize[fields[0]] = int(fields[1])
	return chromSize

def replace_nan(lst):
	for i,v in enumerate(lst):
		if numpy.isnan(v):
			lst[i]=0
	return lst
def all_nan(lst):
	'''check empty list'''
	for i,v in enumerate(lst):
		if not numpy.isnan(v):
			return False
	return True
	
def local_maximum(lst):
	'''pick out local maximum interval'''	#[Interval(119428445, 119428446, value=[3, 37.2617375672]), Interval(119428516, 119428517, value=[8, 84.8253807823])]
	
	mydict={}
	for intval in lst:
		#start_pos = intval.start
		end_pos = intval.end
		#coverage = intval.value[0]
		pvalue = intval.value[1]
		mydict[end_pos] = pvalue
	max_coord = max(mydict, key=mydict.get)
	return(max_coord, mydict[max_coord])		#return coordinate and pvalue
			
def main():
	usage="%prog [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-p","--peak-file",action="store",type="string",dest="peak_file",help="Peak file generated by ChEAP_PeakCalling")
	parser.add_option("-f","--forward",action="store",type="string",dest="forward_peak",help="BigWig file of forward peak (first 5nt)")
	parser.add_option("-r","--reverse",action="store",type="string",dest="reverse_peak",help="BigWig file of reverse peak (first 5nt)")
	parser.add_option("-c","--chromSize",action="store",type="string",dest="chromSize",help="Chromosome size file. Tab or space separated text file with 2 columns: first column is chromosome name, second column is size of the chromosome.")
	parser.add_option("-w","--window",action="store",type="int",dest="window_size",default=5,help="Window size (on genome) to calculate cross strand distance. default=%default")
	parser.add_option("-s","--shift-size",action="store",type="int",dest="max_distance",default=100,help="Maximum shift size. default=%default")

	(options,args)=parser.parse_args()

	if not (options.peak_file and options.forward_peak and options.reverse_peak and options.chromSize):
		parser.print_help()
		sys.exit(0)
	if options.window_size <1:
		print >>sys.stderr, "window size must be intreger larger than 1"
		parser.print_help()
		sys.exit(0)	
	fwd = BigWigFile( file=open(options.forward_peak) )
	rev = BigWigFile( file=open(options.reverse_peak) )
	chrom_sizes = load_chromsize(options.chromSize)	
	shiftSize=collections.defaultdict(list)
	count=0
	avg_eud=collections.defaultdict(int)	#average euclidean distance over window
	for line in open(options.peak_file,'r'):
		if line.startswith('#'):
			continue
		if not line.rstrip():
			continue
		fields=line.rstrip().split()
		if fields[3] == '-':
			continue
		if int(fields[4]) <30:
			continue
		chrom = fields[0]
		peak_pos = int(fields[2])
		peak_start = peak_pos - options.window_size
		peak_end = peak_pos + options.window_size
		if peak_start <0:
			peak_start=0
		if peak_end > chrom_sizes[chrom]:
			peak_end = chrom_sizes[chrom]
		
		fwd_signal = fwd.get_as_array(chrom,peak_start,peak_end)
		#if all_nan(fwd_signal):
		#	continue
		fwd_signal = replace_nan( fwd_signal )
		for offset in range(0,options.max_distance+1):
			rev_signal = rev.get_as_array(chrom,peak_start + offset, peak_end + offset)
			rev_signal = replace_nan( rev_signal )
			#print >>OUT, chrom + ":" + str(peak_start) + '-' + str(peak_end) + '\t' + str(offset) + '\t' + str(twoList.euclidean_distance(fwd_signal,rev_signal))
			shiftSize[chrom + str(peak_pos)].append(twoList.euclidean_distance(fwd_signal,rev_signal))
	for k in shiftSize:
		if len(set(shiftSize[k]))==1:
			continue
		count +=1
		norm_factor = max(shiftSize[k])
		for indx, val in enumerate(shiftSize[k]):
			avg_eud[indx] += val/norm_factor
	for k,v in avg_eud.items():
		print str(k) + '\t' + str(v/count)
			
if __name__=='__main__':
	main()