#!/usr/bin/env python

x = open('knownGene.hg19.bed','r')
prom = open('knownGene_promoter.hg19.bed','w')
prox = open('knownGene_proximal.hg19.bed','w')
dist = open('knownGene_distal.hg19.bed','w')

def promoter_coord(st):
	l = st - 500
	if l < 0:
		l = 1
	return [str(l),str(st+500)]

def proximal_coord(st,side):
	if side == 'left':
		l = st - 10000
		r = st - 500
		if l < 0:
			l = 1
		if r < 0:
			r = 1
		return [str(l),str(r)]
	else:
		return [str(st+500),str(st+10000)]
def distal_coord(st,side):
	if side == 'left':
		l = st - 50000
		r = st - 10000
		if l < 0:
			l = 1
		if r < 0:
			r = 1
		return [str(l),str(r)]
	else:
		return [str(st+10000),str(st+50000)]

for line in x:
	s = line.strip().split('\t')
	if s[5] == '+':
		promo = promoter_coord(int(s[1]))
		lproxi = proximal_coord(int(s[1]),'left')
		rproxi = proximal_coord(int(s[1]),'right')
		ldist = distal_coord(int(s[1]),'left')
		rdist = distal_coord(int(s[1]),'right')
		y = [s[0],promo[0],promo[1],s[3],s[4],s[5],'\n']
		z1 = [s[0],lproxi[0],lproxi[1],s[3],s[4],s[5],'\n']
		z2 = [s[0],rproxi[0],rproxi[1],s[3],s[4],s[5],'\n']
		x1 = [s[0],ldist[0],ldist[1],s[3],s[4],s[5],'\n']
		x2 = [s[0],rdist[0],rdist[1],s[3],s[4],s[5],'\n']
		prom.write('\t'.join(y))
		prox.write('\t'.join(z1))
		prox.write('\t'.join(z2))
		dist.write('\t'.join(x1))
		dist.write('\t'.join(x2))
	else:
		promo = promoter_coord(int(s[2]))
		lproxi = proximal_coord(int(s[2]),'left')
		rproxi = proximal_coord(int(s[2]),'right')
		ldist = distal_coord(int(s[2]),'left')
		rdist = distal_coord(int(s[2]),'right')
		y = [s[0],promo[0],promo[1],s[3],s[4],s[5],'\n']
		z1 = [s[0],lproxi[0],lproxi[1],s[3],s[4],s[5],'\n']
		z2 = [s[0],rproxi[0],rproxi[1],s[3],s[4],s[5],'\n']
		x1 = [s[0],ldist[0],ldist[1],s[3],s[4],s[5],'\n']
		x2 = [s[0],rdist[0],rdist[1],s[3],s[4],s[5],'\n']
		prom.write('\t'.join(y))
		prox.write('\t'.join(z1))
		prox.write('\t'.join(z2))
		dist.write('\t'.join(x1))
		dist.write('\t'.join(x2))

x.close()
prom.close()
prox.close()
dist.close()