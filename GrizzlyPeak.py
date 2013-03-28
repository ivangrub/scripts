#!/usr/bin/env python

import numpy as np

def PeakShape():
	KK = np.floor(int(args.d)/int(args.w))
	mu = int(args.d)
	r = 2.
	p = r/(mu+r)
	X = np.random.negative_binomial(r,p,1e6)
	X2 = np.ceil(X*np.random.rand(1,1e6))
	I  = np.nonzero(np.logical_and(X > 60, X < 170))[1]
	L1 = np.ceil((1-X2[0][I]/args.w)
	L2 = np.ceil((X[I]-X2[0][I])/args.w)
	H1 = np.histogram(L1,np.arange((-1000/int(args.w),0)))
	H2 = np.hist( L2,np.arange(0,(1000/int(args.w))))
	H1 = H1/sum(H1)
	H2 = H2/sum(H2)
	dat1 = [H1 np.zeros(1000/int(args.w))]
	c1 = np.convolve(dat1,np.ones(KK))
	c1 = c1[:-KK+1]
	dat2 = [np.zeros(1000/int(args.w)) H2]
	c2 = np.convolve(dat2,np.ones(KK))
	c2 = c2[KK:]
	c = c1 + c2
	c = smooth(c,5)
	c = c/max(c)
	return c

def Extrap():
	shape = shape[shape > 0]
	lshape = len(shape)
	shp = shape[shape > 0]
	shp = shp/max(shp)
	lshp = len(shp)
	shape2 = shape[shape > 0]
	shape2 = shape2/max(shape2)
	lshape2 = len(shape2)
	shape2 = np.interp(np.arange(1:int(args.w)*lshape2:int(args.w)),shape2,np.arange(1,jump*lshape2))
	lshape2 = len(shape2)

def InitData():

def PFR():

def FitShape():
shape = PeakShape()
Extra(shape)
InitData()