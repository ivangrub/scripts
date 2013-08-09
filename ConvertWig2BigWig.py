#!/usr/bin/env python

import os

names = ['chip.CSEM_ES_Pol2.k100.wig','chip.ES_Pol2_150bp_neigh1_wB.mm9_norandom.wig','chip.ES_Pol2.k100.mm9_norandom.wig','chip.TY_ES_PolII.mm9_norandom.m1.mm9_norandom.wig']

for type in names:
	new = type.split('.wig')
	os.system('wigToBigWig %s mm9_norandom_size.txt %s.bigwig' % (type,new[0]))