#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
Convert BAM/SAM file into wig file. BAM file must be sorted and indexed
How to sort and index BAM file: http://genome.ucsc.edu/goldenPath/help/bam.html
-------------------------------------------------------------------------------------------------'''

#import built-in modules
import os,sys
import string
from optparse import OptionParser
import warnings
import string
import collections

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

#import my own modules
from qcmodule import SAM
from qcmodule import BED
#changes to the paths

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyright 2012, Wei Li's Lab"
__credits__ = []
__license__ = "GPL"
__version__="2.0"
__maintainer__ = "Liguo Wang"
__email__ = "liguow@bcm.edu"
__status__ = "Development" #Prototype or Production


def printlog (mesg):
	'''print progress into stderr and log file'''
	mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	LOG=open('class.log','a')
	print >>sys.stderr,mesg
	print >>LOG,mesg

def load_chromsize(file):
	'''read chrom.size file'''
	chromSize={}
	for line in open(file,'r'):
		if line.startswith('#'):continue
		if not line.strip():continue
		fields = line.strip().split()
		chromSize[fields[0]] = int(fields[1])
	return chromSize

def build_wig(obj,ext):
	'''build wig files from alignedReads object'''
	
	Pwig=collections.defaultdict(int)
	Nwig=collections.defaultdict(int)
	if obj is None:
		return (Pwig,Nwig)
	for aligned_read in obj:
		flag=0
		if aligned_read.is_qcfail:			#skip low quanlity
			continue
		if aligned_read.is_duplicate:		#skip duplicate read
			continue
		if aligned_read.is_secondary:		#skip non primary hit
			continue
		if aligned_read.is_unmapped:		#skip unmap read
			continue		
		read_tags = aligned_read.tags
		for i in read_tags:
			if i[0] in SAM.ParseBAM.multi_hit_tags and i[1] >1:
				flag=1						#multiple hit read
				break
		if flag==1:
			continue					#skip multiple map read				
		
		if ext is None: ext = aligned_read.qlen
		if aligned_read.is_reverse:
			read5end = aligned_read.pos + aligned_read.qlen
			read3end = read5end - ext
			for i in xrange(read3end +1, read5end +1):
				Nwig[i] += 1
		else:
			read5end = aligned_read.pos
			read3end = read5end + ext
			for i in xrange(read5end +1, read3end +1):
				Pwig[i] += 1
	return (Pwig,Nwig)
			
			
def main():
	usage="%prog [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Input file in BAM format. BAM file must be sorted and indexed using samTools. HowTo: http://genome.ucsc.edu/goldenPath/help/bam.html")
	parser.add_option("-r","--chromSize",action="store",type="string",dest="chromSize",help="Chromosome size file. Tab or space separated text file with 2 columns: first column is chromosome name, second column is size of the chromosome.")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output wig files(s). \"Prefix_Forward.wig\" and \"Prefix_Reverse.wig\" will be generated")
	parser.add_option("-b","--bin",action="store",type="int",dest="bin",default=100000,help="Chromosome chunk size. Each chomosome will be cut into samll chunks of this size. Decrease chunk size will save more RAM. default=%default (bp)")
	parser.add_option("-e","--extension",action="store",type="int",dest="extension",default=None,help="Extended coverage from 5' end of read. default=%default (full read coverage will be used)")

	(options,args)=parser.parse_args()

	if not (options.output_prefix and options.input_file and options.chromSize):
		parser.print_help()
		sys.exit(0)
	for file in (options.input_file,options.chromSize):
		if not os.path.exists(file):
			print >>sys.stderr, '\n\n' + file + " does NOT exists" + '\n'
			sys.exit(0)
	if not os.path.exists(options.input_file + '.bai'):
		print >>sys.stderr, "index file " + options.input_file + '.bai' + "does not exists"
		sys.exit(0)


	chrom_sizes = load_chromsize(options.chromSize)
	samfile = SAM.ParseBAM(options.input_file)
	FWOUT = open(options.output_prefix + "_Forward.wig",'w')
	RWOUT = open(options.output_prefix + "_Reverse.wig",'w')
	
	for chr_name, chr_size in chrom_sizes.items():		#iterate each chrom
		try:
			samfile.fetchAlignments(chr_name,0,chr_size)
		except:
			print >>sys.stderr, "No alignments for " + chr_name + '. skipped'
			continue
		print >>sys.stderr, "Processing " + chr_name + " ..."
		FWOUT.write('variableStep chrom='+chr_name+'\n')
		RWOUT.write('variableStep chrom='+chr_name+'\n')
		for interval in BED.tillingBed(chrName = chr_name,chrSize = chr_size,stepSize = options.bin):	#cut chrom into bins, interval such as ('chr1', 235000000, 236000000)
			Fwig={}
			Rwig={}
			alignedReads = samfile.fetchAlignments(interval[0],interval[1],interval[2])
			(Fwig,Rwig) = build_wig(alignedReads,options.extension)
			
			if (len(Fwig)>0):
				for i in xrange(interval[1]+1,interval[2]+1):
					if Fwig.has_key(i):
						FWOUT.write("%d\t%d\n" % (i, Fwig[i])) 
			if (len(Rwig)>0):
				for i in xrange(interval[1]+1,interval[2]+1):
					if Rwig.has_key(i):
						RWOUT.write("%d\t%d\n" % (i, Rwig[i]))

if __name__ == '__main__':
	main()
