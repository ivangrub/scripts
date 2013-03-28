#!/usr/bin/env python

import argparse
import pysam as pys
from Bio import SeqIO
import os

parser = argparse.ArgumentParser(description='Process an annotated BAM file into global coordinates')
parser.add_argument('-i', type=str,help='Input BAM file')
parser.add_argument('-f',type = str, help='Reference annotation')
parser.add_argument('-s',type = str, help='Chromosome lengths')
parser.add_argument('-o',type=str,help='Output BAM file')

args = parser.parse_args()

dir = os.environ['BOWTIE_INDEXES']
fasta = open(dir+'/'+args.f,'r')
size = open(dir+'/'+args.s,'r')

ref = []
header = {}
header['HD'] = {'VN':'1.0'}
for i in size:
	b = i.strip().split('\t')
	ref.append({'LN':int(b[1]), 'SN':b[0]})
header['SQ'] = ref

chrdir = {}
for seq in SeqIO.parse(fasta,'fasta'):
	info = seq.description
	b = info.split(' ')
	chrdir[b[0]] = b[1]	

x = pys.Samfile(args.i,'rb')
out = pys.Samfile('tmp.sam','wh',header = header)
bamout = pys.Samfile(args.o,'wb',template = out)

ref = bamout.references
reflist = {}
j = 0
for i in ref:
	reflist[i] = j
	j += 1

for seq in x:
	if seq.is_unmapped:
		continue
	name = chrdir[str(seq.rname)]
	bin = name.split(':')
	chr = bin[0]
	st = bin[1].split('-')
	pos_init = int(seq.pos)
	pos_new = pos_init + int(st[0]) - 1
	seq.pos = pos_new
	seq.rname = reflist[chr]
	bamout.write(seq)

x.close()
out.close()
bamout.close()
