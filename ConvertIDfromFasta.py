#!/usr/bin/env python

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Convert target-id back to transcript names.')
parser.add_argument('-ref',help = 'Enter the name of the transcript fasta file. The default is ens_transcripts.hg19.fa',default = 'ens_transcripts.hg19.fa')
parser.add_argument('-o',help = 'Input the results.xprs file',default = 'results.xprs')
args = parser.parse_args()

fasta_name = open(args.ref,'r')
changeid = open(args.o,'r')
modid = open('results_mod.xprs','w')

convert = {}
for seq in SeqIO.parse(fasta_name,'fasta'):
	id_info = seq.description.strip().split(' ')
	convert[id_info[0]] = id_info[1]

fasta_name.close()

i = 0	
for line in changeid:
	if i > 0:
		s = line.strip().split('\t')
		modid.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'  % (s[0],convert[s[1]],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[12],s[13]))
		i += 1
	else:
		modid.write(line)
		i += 1
		
changeid.close()
modid.close()	

