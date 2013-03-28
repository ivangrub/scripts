#!/usr/bin/env python

import os

dir = '/Users/ivang/Desktop/CC_Rad23b_PostSPP'
for i in os.listdir('%s' % dir):
	print 'Adding %s' % i
	os.system('cat %s/%s >> %s/Rad23b_PostSPP.wig' % (dir,i,dir))
	