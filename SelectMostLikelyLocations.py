#!/usr/bin/env python	

import pysam as pys
import argparse

parser = argparse.ArgumentParser(description='Project SAM/BAM format reads back onto the appropriate genome and call significant peaks.')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = '')
parser.add_argument('-o',help ='Name of outputfile',default = '')
args = parser.parse_args()

r = pys.Samfile(args.r,'rb')
out = pys.Samfile('Mostlikely_%s.bam' % args.o,'wb',template = r)

nreads = 0
odds = []
reads = []

for seq in r:
	if nreads == 0:
		nreads += 1
		prev_seq = seq
		frac_str = seq.opt('XP')
		frac = float(frac_str)
		odds.append(frac)
		reads.append(seq)
		
		continue
		
	if prev_seq.qname == seq.qname:
		prev_seq = seq
		frac_str = seq.opt('XP')
		frac = float(frac_str)
		odds.append(frac)
		reads.append(seq)
		
		continue
	else:
		ind = odds.index(max(odds))
		
		out.write(reads[ind])
		odds = []
		reads = []
		
		prev_seq = seq
		frac_str = seq.opt('XP')
		frac = float(frac_str)
		odds.append(frac)
		reads.append(seq)
		
out.close()
r.close()
		
	