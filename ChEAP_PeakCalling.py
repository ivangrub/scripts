#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
CHEAP: ChIP-exo Analysis Program
-------------------------------------------------------------------------------------------------'''

#import built-in modules
import os,sys
import string
from optparse import OptionParser
import warnings
import string
import collections
import  numpy
from time import strftime
import signal
import math


#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from bx.bbi.bigwig_file import BigWigFile

#import my own modules
from qcmodule import SAM
from qcmodule import BED
from qcmodule import poisson

#changes to the paths

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyright 2012, Baylor College of Medicine"
__credits__ = []
__license__ = "GPL"
__version__="2.0"
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
ChIP-Exo  Analysis  Pipeline (peak calling)                                      
'''                                                                                  

                                                                                        
def printlog (mesg):
	'''print progress into stderr and log file'''
	mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	#LOG=open('ChEAP_running.log','a')
	print >>sys.stderr,mesg
	#print >>LOG,mesg

def signal_handler(signal, frame):
	print >>sys.stderr, '\nYou pressed Ctrl+C. ChEAP was terminated!'
	sys.exit(0)

def load_chromsize(file):
	'''read chrom.size file'''
	chromSize={}
	for line in open(file,'r'):
		if line.startswith('#'):continue
		if not line.strip():continue
		fields = line.strip().split()
		chromSize[fields[0]] = int(fields[1])
	return chromSize

def peakRoot(obj):
	'''Return coverage for each peak root. cordindate is 1 based'''
	
	peak_root=collections.defaultdict(int)
	if obj is None:
		return peak_root
	for aligned_read in obj:
		
		#Input BAM file is supposed to have good quality, unique reads
		#flag=0
		#if aligned_read.is_qcfail:			#skip low quanlity
		#	continue
		#if aligned_read.is_duplicate:		#skip duplicate read
		#	continue
		#if aligned_read.is_secondary:		#skip non primary hit
		#	continue
		#if aligned_read.is_unmapped:		#skip unmap read
		#	continue		
		#read_tags = aligned_read.tags
		#for i in read_tags:
		#	if i[0] in SAM.ParseBAM.multi_hit_tags and i[1] >1:
		#		flag=1						#multiple hit read
		#		break
		#if flag==1:
		#	continue					#skip multiple map read				
		
		if aligned_read.is_reverse:
			strand = '-'
			read5end = aligned_read.pos + aligned_read.qlen
		else:
			strand = '+'
			read5end = aligned_read.pos + 1
		peak_root[str(read5end) + ":" + strand] +=1
	return peak_root

def peak_area(peak_root, chrom, ext, fw_bw,re_bw, fold_cut):
	'''peak_root is dict returned by peakRoot, ext=number of nucleotide extended from 5' end of read
	fw_bw =BigWigFile object of forward strand rv_bw = BigWigFIle object of reverse strand, fold_cut
	= coverage cutoff of a valid peak root'''
	
	peak_area = collections.defaultdict(int)
	for key in peak_root:
		if peak_root[key] < fold_cut:	#skip if this peak root does NOT meet our criteria
			continue
		(peak_pos, peak_strand) = key.split(':')
		
		if peak_strand == '-':
			region_end = int(peak_pos)
			region_start = region_end - ext
			if region_start < 0:region_start = 0
			obj = re_bw
		if peak_strand == '+':
			region_start = int(peak_pos)-1
			if region_start < 0:region_start = 0
			region_end = region_start + ext
			obj = fw_bw
				
		for val in obj.get_as_array(chrom,region_st,region_end):
			if numpy.isnan(val):continue
			peak_area[key] += val
	return peak_area
			

def sum_bwfile(chrom,root_pos,ext_size,bw_obj,chrmSize):
	'''return maximum sum signal for a specified region. root_pos is 1 based'''
	
	sum_list=[]
	for start in range(root_pos - ext_size, root_pos):
		st = start
		if st<0: st=0
		end = st + ext_size
		if end > chrmSize[chrom]: end = chrmSize[chrom]
		wig_sum = numpy.nansum(bw_obj.get_as_array(chrom,st,end))
		sum_list.append(wig_sum)
	return max(sum_list)

def merge_peaks(input_dict,fuzziness=5):
	'''merge closest peaks together'''
	
	d={}
	for k in input_dict.keys():
		flag=0
		(chrom,coord,strand) = k.split("\t")
		coord = int(coord)
		for i in range(coord - fuzziness,coord + fuzziness +1):
			k2 = chrom + '\t' + str(i) + '\t' +strand
			if k2 == k:
				continue
			if (k2 in input_dict) and (input_dict[k2] > input_dict[k]):
				flag=1
				break
		if flag==0:
			d[k] = input_dict[k]
	return d

def cal_poisson_pvalue(q,start,end, tree_obj,window_size,bg_root_num):
	'''calculate poisson pvalue'''
	
	#pos2val={}
	bg_signal=[]
	
	bg_start = end - int(window_size/2) - 1		#backkground window start
	bg_end = end + int(window_size/2)				#background window end
	all_intervals = tree_obj.find(bg_start,bg_end)
	
	if (len(all_intervals) >= bg_root_num):
		for intv in all_intervals:
			bg_signal.append(intv.value)
			lamda = numpy.median(bg_signal)			#this is local mean
	else:
		up_intervals = tree_obj.upstream_of_interval(Interval(start,end,strand='+'),num_intervals=bg_root_num)
		down_intervals = tree_obj.downstream_of_interval(Interval(start,end,strand='+'),num_intervals=bg_root_num)
		for intv in up_intervals:bg_signal.append(intv.value)	
		for intv in down_intervals:bg_signal.append(intv.value)
		lamda = numpy.median(bg_signal)			#this is local mean
	if len(bg_signal)==0:
		return 0
	if lamda > q:
		return 0
	return poisson.cumu_poip(q, lamda, logp=True)
		
def main():
	usage="%prog [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-b","--forward",action="store",type="string",dest="forward_bw",help="BigWig file for forward reads (extend 1 nt from 5' end of read)")
	parser.add_option("-d","--reverse",action="store",type="string",dest="reverse_bw",help="BigWig file for reverse reads (extend 1 nt from 5' end of read)")
	parser.add_option("-s","--chromSize",action="store",type="string",dest="chromSize",help="Chromosome size file. Tab or space separated text file with 2 columns: first column is chromosome name, second column is size of the chromosome.")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files")
	parser.add_option("-z","--fuzziness",action="store",type="int",dest="fuzzy_size",default=10,help="Peaks within fuzzy window will be merged. default=%default (bp)")
	parser.add_option("-w","--bgw",action="store",type="int",dest="window_size",default=200,help="Background window size used to determine background signal level (lambda in Poisson model). default=%default (bp)")
	parser.add_option("-c","--chunk",action="store",type="int",dest="chunk_size",default=100000,help="Chromosome chunk size. Each chomosome will be cut into samll chunks of this size. Decrease chunk size will save more RAM. default=%default (bp)")
	parser.add_option("-p","--pvalue",action="store",type="float",dest="pvalue_cutoff",default=0.1,help="Pvalue cutoff for peak detection. default=%default")
	parser.add_option("-r","--bg-root-num",action="store",type="float",dest="bg_root_num",default=100,help="Background peak root number. default=%default")
	parser.add_option("-e","--extention",action="store",type="int",dest="extention_size",default=5,help="Window size used to calculate peak area. Larger number will signficantly reduce speed, and make peak calling more meaningless.  default=%default")

	(options,args)=parser.parse_args()

	if not (options.output_prefix and options.chromSize and options.forward_bw and options.reverse_bw):
		parser.print_help()
		sys.exit(0)
	for file in (options.chromSize,options.forward_bw,options.reverse_bw):
		if not os.path.exists(file):
			print >>sys.stderr, '\n\n' + file + " does NOT exists" + '\n'
			sys.exit(0)
	
	chrom_sizes = load_chromsize(options.chromSize)
	OUT = open(options.output_prefix + ".single_nt_peak.xls",'w')
	fw_bw_obj = BigWigFile( file = open(options.forward_bw))
	rv_bw_obj = BigWigFile( file = open(options.reverse_bw))
	rv_peak_roots = {}
	rv_peak_height = {}
	rv_ranges={}
	rv_peak_pvalue={}
	pv_cutoff = -10*math.log10(options.pvalue_cutoff)	
	signal.signal(signal.SIGINT, signal_handler)


	print >>sys.stderr, logo	
	
	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#calculate peak height and peak area for forward bigwig
	print >>sys.stderr, "@" + strftime("%Y-%m-%d %H:%M:%S") + ": Processing " + options.forward_bw + '  ...'
	for chr_name, chr_size in chrom_sizes.items():		#iterate each chrom
		fw_peak_roots = {}	#key is chr,pos,strand,height: ("chr19   51345387        +       2.83"), value is area("2.82999992371")
		fw_peak_height = {}
		fw_ranges={}
		fw_peak_pvalue={}
		if chr_name != 'chrY':
			continue
		print >>sys.stderr, "@" + strftime("%Y-%m-%d %H:%M:%S") + ": Processing " + chr_name + " ..."
		progress = 0
		coord = 0	
		#for each chunk
		for interval in BED.tillingBed(chrName = chr_name,chrSize = chr_size,stepSize = options.chunk_size):	#cut chrom into bins, interval such as ('chr1', 235000000, 236000000)				
			for indx,val in enumerate(fw_bw_obj.get_as_array(interval[0],interval[1],interval[2])):
				coord += 1	#coord is 1-based on genome
				if numpy.isnan(val):continue
				area_value = sum_bwfile(chr_name, coord, options.extention_size, fw_bw_obj,chrom_sizes)
				fw_peak_roots[chr_name + "\t" + str(coord) + "\t+"] = area_value		#key is chrom + position + strand,value is area
				fw_peak_height[chr_name + "\t" + str(coord) + "\t+"] = val
				if chr_name not in fw_ranges:
					fw_ranges[chr_name] = IntervalTree()
				else:
					fw_ranges[chr_name].insert_interval( Interval( coord-1, coord, value=area_value) )
			finish_part = int(interval[2]*100/chr_size)
			if finish_part > progress:
				print >>sys.stderr, " %d%% finished\r" % (finish_part),
				progress = finish_part	
	
	
		#fw_global_lamda = numpy.mean(fw_peak_roots.values())
		#print >>sys.stderr, "Global mean (Forward) = " + str(fw_global_lamda)
		print >>sys.stderr, "@" + strftime("%Y-%m-%d %H:%M:%S") + ": Calculating pvalues for " + options.forward_bw + '  ...'
		for k in fw_peak_roots:
			chrom = k.split("\t")[0]
			coord = int(k.split("\t")[1])
			fw_peak_pvalue[k] = cal_poisson_pvalue(int(fw_peak_roots[k]), coord-1, coord, fw_ranges[chrom],options.window_size,options.bg_root_num)

	
		fw_peak_filtered = merge_peaks(fw_peak_height,fuzziness=options.fuzzy_size)	
		for k,v in fw_peak_filtered.items():
			#print k + '\t' + str(v)
			(chrom,end,strand) = k.split('\t')
			end = int(end)
			start = end -1
			height = str(v)
			area = str(fw_peak_roots[k])
			pvalue = fw_peak_pvalue[k]
			if pvalue < pv_cutoff:continue
			print >>OUT, '\t'.join([chrom, str(start), str(end), area,str(round(pvalue)),strand,height])
	
	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#calculate peak height and peak area for reverse bigwig
	print >>sys.stderr, "@" + strftime("%Y-%m-%d %H:%M:%S") + ": Processing " + options.reverse_bw + '  ...'
	for chr_name, chr_size in chrom_sizes.items():		#iterate each chrom
		if chr_name != 'chrY':
			continue
		print >>sys.stderr, "@" + strftime("%Y-%m-%d %H:%M:%S") + ": Processing " + chr_name + " ..."
		progress = 0
		coord = 0	
		#for each chunk
		for interval in BED.tillingBed(chrName = chr_name,chrSize = chr_size,stepSize = options.chunk_size):	#cut chrom into bins, interval such as ('chr1', 235000000, 236000000)				
			
			for indx,val in enumerate(rv_bw_obj.get_as_array(interval[0],interval[1],interval[2])):
				coord += 1	#coord is 1-based on genome
				if numpy.isnan(val):continue
				area_value = sum_bwfile(chr_name, coord, options.extention_size, rv_bw_obj,chrom_sizes)				
				rv_peak_roots[chr_name + "\t" + str(coord) + "\t-"] = area_value
				rv_peak_height[chr_name + "\t" + str(coord) + "\t-"] = val
				if chr_name not in rv_ranges:
					rv_ranges[chr_name] = IntervalTree()
				else:
					rv_ranges[chr_name].insert_interval( Interval( coord-1, coord, value = area_value) )
			finish_part = int(interval[2]*100/chr_size)
			if finish_part > progress:
				print >>sys.stderr, " %d%% finished\r" % (finish_part),
				progress = finish_part
	

	#rv_global_lamda = numpy.mean(rv_peak_roots.values())
	#print >>sys.stderr, "Global mean (Reverse) = " + str(rv_global_lamda)
	print >>sys.stderr, "@" + strftime("%Y-%m-%d %H:%M:%S") + ": Calculating pvalues for " + options.reverse_bw + '  ... '
	for k in rv_peak_roots:
		chrom = k.split("\t")[0]
		coord = int(k.split("\t")[1])
		rv_peak_pvalue[k] = cal_poisson_pvalue(int(rv_peak_roots[k]),coord-1,coord, rv_ranges[chrom],options.window_size,options.bg_root_num)
		#print k + '\t' + str(rv_peak_roots[k]) + '\t' + str(pvalue)


	rv_peak_filtered = merge_peaks(rv_peak_height,fuzziness=options.fuzzy_size)
	for k,v in rv_peak_filtered.items():
		(chrom,end,strand) = k.split('\t')
		end = int(end)
		start = end -1
		height = str(v)
		area = str(rv_peak_roots[k])
		pvalue = rv_peak_pvalue[k]
		if pvalue < pv_cutoff:continue
		
		print >>OUT, '\t'.join([chrom, str(start), str(end), area, str(round(pvalue)),strand,height])


if __name__ == '__main__':
	main()
