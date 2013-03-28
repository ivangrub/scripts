#!/usr/bin/env python

import numpy as np 
import argparse

def readpeaks(file):
	o = open(file,'r')
	peak = []
	for line in o:
		s = line.strip().split()
		peak.append([s[0],s[1],s[2]])
	return np.array(peak)

def left(sub,ref):

def right(sub,ref):

def inside(sub,ref):

def contain(sub,ref):


parser = argparse.ArgumentParser(description='Establish overlap of any peak calling method for ChIP-Seq')
parser.add_argument('-p1',help ='Read subject peak file.',default = '')
parser.add_argument('-p2',help =  'Read reference peak file.',default = '')
parser.add_argument('-m',help = 'Minimum basepair overlap',default = 1)
parser.add_argument('-o',help = 'The prefix of the output files.',default = 'test')
args = parser.parse_args()

subject = readpeaks(args.p1)
reference = readpeaks(args.p2)

overlapped = open('%s_Overlapped.txt' % args.o,'w')
unique = open('%s_Unique.txt' % args.o,'w')

for p in subject:
	if left(p,reference) or right(p,reference) or inside(p,reference) or contain(p,reference):
		overlapped.write('%s\t%d\t%d\n' % (p[0],p[1],p[2]))
	else:
		unique.write('%s\t%d\t%d\n' % (p[0],p[1],p[2]))




