# Initialize the classes and make the Peak/Overlap objects
#
# By: Ivan Grubisic
# Date: December 6th, 2011
#

from numpy import *
from matplotlib.pyplot import *
from scipy.stats import *

# Define the peak class
class Peak(object):
	"""Defines the peak object which saves the chromosome, left and right edge coordinates of 
	the peak, the peak's p-value and calculates the midpoint and length of the peak"""
	
	def __init__(self,file):
		self.chr = []						# The chromosome that the peak is on
		self.left = []						# The left edge of the peak
		self.right = []						# The right edge of the peak
		self.enrich = []					# The IP enrichment over background
		self.mid = []						# The midpoint of the peak
		self.peaklen = []					# The length of the peak

		# Iterate through the text file
		for line in file:
			s = line.strip().split()
			self.chr.append(s[0])
			self.left.append(float(s[1]))
			self.right.append(float(s[2]))
			self.enrich.append(float(s[6]))
			self.mid.append((float(s[1]) + float(s[2]))/2.)
			self.peaklen.append(float(s[2]) - float(s[1]))

	def hist_peaklen(self,iplist):
		"""Create a histogram plot of the peak length"""
		hist(self.peaklen, bins = 100)
		xlabel('Peak Length (bp)')
		ylabel('Count')
		
		name = iplist[0]
		for j in range(len(iplist)-1):
			name = name + '_' + iplist[j+1]
		title(name)
		grid(False)
		show()
			
	def upstream(self,right,left,index):
		"""Determine if reference peak is on the left edge of the second IP"""
		ind = (self.left[index] < left) & (self.right[index] > left) & (self.right[index] <= right)
		return ind
		
	def inside(self,right,left,index):
		"""Determine if the reference peak fits inside the second IP"""
		ind = (self.left[index] >= left) & (self.right[index] <= right) 
		return ind
		
	def other_inside(self,right,left,index):
		"""Determine if the second IP fits inside the reference peak"""
		ind = (self.left[index] <= left) & (self.right[index] > right)
		return ind
		
	def downstream(self,right,left,index):
		"""Determine if the reference peak is on the right edge of the second IP"""
		ind = (self.left[index] <= right) & (self.right[index] > right) & (self.left[index] >= left)
		return ind

# Define the overlap class and the ability for multiple overlaps	
class Overlap(Peak):
	"""Create a class object of overlaps between 2 IPs"""
	def __init__(self,p1,p2,percent = 1):
		print 'Initializing the Overlap and save the information regarding the non-overlapping peaks'
		self.chr = []
		self.left = []
		self.right = []
		self.mid = []
		self.enrich1 = []
		self.enrich2 = []
		self.peaklen = []
		self.interdist = []
		self.nochr = []
		self.noleft = []
		self.noright = []
		self.nomid = []
		self.noenrich1 = []
		self.nopeaklen = []
		
		# Take the top X percent of peaks to calculate the overlap
		num1 = int(len(p1.chr)*percent)
		num2 = int(len(p2.chr)*percent)
		
		c = array(p2.chr[:num2])
		right = array(p2.right[:num2])
		left = array(p2.left[:num2])
		enrich = array(p2.enrich[:num2])
		mid = array(p2.mid[:num2])
		
		for j in range(len(p1.left[:num1])):
			# Look at only the peaks that are on the same chromosome
			chr = p1.chr[j] == c
			r = right[chr]
			l = left[chr]
			midpoint = mid[chr]
			p = enrich[chr]
			
			# Determine if there is any overlap at all
			ind = p1.upstream(r,l,j) | p1.inside(r,l,j) | p1.other_inside(r,l,j) | p1.downstream(r,l,j)
			if ind.any() == False:
				self.nochr.append(p1.chr[j])
				self.noleft.append(p1.left[j])
				self.noright.append(p1.right[j])
				self.nomid.append(p1.mid[j])
				self.noenrich1.append(p1.enrich1[j])
				self.nopeaklen.append(p1.peaklen[j])
			
			# Find the index of the overlapping peaks within the refernce
			x = where(ind)[0]
			
			# Create a new object with the same qualities for all of the possible overlaps
			for i in x:
				self.chr.append(p1.chr[j])
				self.left.append(min(p1.left[j],l[i]))
				self.right.append(max(p1.right[j],r[i]))
				self.mid.append(mean([self.left,self.right]))
				self.enrich1.append(p1.enrich[j])
				self.enrich2.append(p[i])
				self.peaklen.append(self.right[-1] - self.left[-1])
			
			# Return the average distance between the peak midpoints. 
			# The length of this list will be equivalent to the number of reference
			# peaks that have an overlap.
			self.interdist.append(mean(p1.mid[j] - midpoint[ind]))
		
	def cross_corr(self,ip):
		"""Plot the enrichment values of the corresponding peak overlaps"""
		
		# Calculate the Pearson coefficient using the scipy module
		r, p = pearsonr(self.enrich1,self.enrich2)
		loglog(self.enrich1,self.enrich2,'b.',lw=0)
		xlabel('%s Enrichment' % ip[0])
		ylabel('%s Enrichment' % ip[-1])
		title('A Pearson Correlation Coefficient of %.2f' % r)
		show()
		
	def fraction(self,p1,percent = 1):
		"""Calculates the number and frac1tion of peaks relative to the reference IP"""
		len_p1 = len(p1.chr)*percent
		len_self = len(self.interdist)
		return len_self, len_self/float(len_p1)*100
	
	def dist_plot(self,iplist):
		"""Plot the distribution of distance between the midpoints of overlapping peaks"""
		n, bins, patches = hist(self.interdist,bins = 100)
		
		# Calculate the number of peaks that fit within 200bp
		num = sum(n[abs(bins[1:]) <= 200])
		
		# Calculate the percentage of possible peaks that fit within 200bp of one another
		per = sum(n[abs(bins[1:]) <= 200])/float(len(self.interdist))*100.
		
		xlabel('Distance (bp)')
		ylabel('Number of Peaks')
		name = iplist[0]
		for j in range(len(iplist)-1):
			name = name + '_' + iplist[j+1]
		title('%s- Num of Peaks within 200bp = %d (%f percent of total)' % (name,num,per))
		show()
		return num, per	
			
class Overlap_NoMerge(Overlap):
	"""Create a class object of overlaps between 2 IPs without merging the windows"""
	def __init__(self,p1,p2,percent = 1):
		print 'Initializing the Overlap'
		self.chr = []
		self.left = []
		self.right = []
		self.mid = []
		self.enrich1 = []
		self.enrich2 = []
		self.peaklen = []
		self.interdist = []
		
		# Take the top X percent of peaks to calculate the overlap
		num1 = int(len(p1.chr)*percent)
		num2 = int(len(p2.chr)*percent)
		
		c = array(p2.chr[:num2])
		right = array(p2.right[:num2])
		left = array(p2.left[:num2])
		enrich = array(p2.enrich[:num2])
		mid = array(p2.mid[:num2])
		
		for j in range(len(p1.left[:num1])):
			# Look at only the peaks that are on the same chromosome
			chr = p1.chr[j] == c
			r = right[chr]
			l = left[chr]
			midpoint = mid[chr]
			p = enrich[chr]
			
			# Determine if there is any overlap at all
			ind = p1.upstream(r,l,j) | p1.inside(r,l,j) | p1.other_inside(r,l,j) | p1.downstream(r,l,j)
			if ind.any() == False:
				continue
			
			# Find the index of the overlapping peaks within the refernce
			x = where(ind)[0]
			
			# Create a new object with the same qualities for all of the possible overlaps
			for i in x:
				self.chr.append(p1.chr[j])
				self.left.append(p1.left[j])
				self.right.append(p1.right[j])
				self.mid.append(mean([self.left,self.right]))
				self.enrich1.append(p1.enrich[j])
				self.enrich2.append(p[i])
				self.peaklen.append(self.right[-1] - self.left[-1])
			
			# Return the average distance between the peak midpoints. 
			# The length of this list will be equivalent to the number of reference
			# peaks that have an overlap.
			self.interdist.append(mean(p1.mid[j] - midpoint[ind]))

class AddOverlap_NoMerge(Overlap):
	"""Create a class object of overlaps between 2 IPs"""
	def __init__(self,p1,p2,percent = 1):
		
		self.chr = []
		self.left = []
		self.right = []
		self.mid = []
		self.enrich1 = []
		self.enrich2 = []
		self.peaklen = []
		self.interdist = []
		
		# Take the top X percent of peaks to calculate the overlap
		num1 = int(len(p1.chr)*percent)
		num2 = int(len(p2.chr)*percent)
		
		c = array(p2.chr[:num2])
		right = array(p2.right[:num2])
		left = array(p2.left[:num2])
		enrich = array(p2.enrich[:num2])
		mid = array(p2.mid[:num2])
		
		for j in range(len(p1.left[:num1])):
			# Only look to overlap peaks on the same chromosome
			chr = p1.chr[j] == c
			r = right[chr]
			l = left[chr]
			midpoint = mid[chr]
			p = enrich[chr]
			
			# Is there any available overlap
			ind = p1.upstream(r,l,j) | p1.inside(r,l,j) | p1.other_inside(r,l,j) | p1.downstream(r,l,j)
			if ind.any() == False:
				continue
			x = where(ind)[0]
			for i in x:
				self.chr.append(p1.chr[j])
				self.left.append(p1.left[j])
				self.right.append(p1.right[j])
				self.mid.append(mean([self.left,self.right]))
				self.enrich1.append(p1.enrich1[j])
				self.enrich2.append(p[i])
				self.peaklen.append(self.right[-1] - self.left[-1])
			self.interdist.append(mean(p1.mid[j] - midpoint[ind]))
			
class AddOverlap(Overlap):
	"""Create a class object of overlaps between 2 IPs"""
	def __init__(self,p1,p2,percent = 1):
		
		self.chr = []
		self.left = []
		self.right = []
		self.mid = []
		self.enrich1 = []
		self.enrich2 = []
		self.peaklen = []
		self.interdist = []
		
		# Take the top X percent of peaks to calculate the overlap
		num1 = int(len(p1.chr)*percent)
		num2 = int(len(p2.chr)*percent)
		
		c = array(p2.chr[:num2])
		right = array(p2.right[:num2])
		left = array(p2.left[:num2])
		enrich = array(p2.enrich[:num2])
		mid = array(p2.mid[:num2])
		
		for j in range(len(p1.left[:num1])):
			# Only look to overlap peaks on the same chromosome
			chr = p1.chr[j] == c
			r = right[chr]
			l = left[chr]
			midpoint = mid[chr]
			p = enrich[chr]
			
			# Is there any available overlap
			ind = p1.upstream(r,l,j) | p1.inside(r,l,j) | p1.other_inside(r,l,j) | p1.downstream(r,l,j)
			if ind.any() == False:
				continue
			x = where(ind)[0]
			for i in x:
				self.chr.append(p1.chr[j])
				self.left.append(min(p1.left[j],l[i]))
				self.right.append(max(p1.right[j],r[i]))
				self.mid.append(mean([self.left,self.right]))
				self.enrich1.append(p1.enrich1[j])
				self.enrich2.append(p[i])
				self.peaklen.append(self.right[-1] - self.left[-1])
			self.interdist.append(mean(p1.mid[j] - midpoint[ind]))			

def annotate_overlap(peak,ip,frac,gene):
	"""Print a text file that gives a gene annotation for the overlapped peaks. It will print 
	the closest TSS and TES, up- and downstream of the middle of the peak window"""
	name = ip[0]
	for j in range(len(ip)-1):
		name = name + '_' + ip[j+1]
	
	print 'Printing the overlapped gene annotation table'
	table = open('Overlap_%s_Top_%sPercent_Peaks.txt' % (name,str(frac)),'w')
	table.write('Gene Name\tGene ID\tChromosome\tStrand\tTSS\tTES\t\tPeak Left Edge\tPeak Right Edge\tReference IP Enrichment\tOverlapping IP Enrichment\n')
	
	# Create numpy arrays for faster computing
	chr = array(peak.chr)
	mid = array(peak.mid)
	
	for k in range(len(peak.chr)):
		# Only look at the genes on the same chromosome as the peak
		c = gene[:][:,4] == chr[k]
		
		# Output the index within the gene array of the appropriate genes
		c_index = c.nonzero()[0]
		
		# Return the TSS and TES values for the available genes
		tss = return_tss(gene[c,:])
		tes = return_tes(gene[c,:])
		
		# Calculate the distance to each TSS/TES from the peak midpoint
		dist_tss = tss - mid[k]
		dist_tes = tes - mid[k]
		
		# Find index of the closest upstream and downstream TSS and TES and return closest gene
		try:
			up_tss = dist_tss[dist_tss >= 0].argmin()
			gene_tss_up = gene[c_index[(dist_tss[dist_tss >= 0][up_tss] == dist_tss).nonzero()[0]]][0]
			table.write('%s\t%s\t%s\t%s\t%s\t%s\t\t%d\t%d\t%f\t%f\n' % (gene_tss_up[3],gene_tss_up[2],gene_tss_up[4],gene_tss_up[5],gene_tss_up[6],gene_tss_up[7],peak.left[k],peak.right[k],peak.enrich1[k],peak.enrich2[k]))
		except ValueError:
			pass
		try:
			down_tss = abs(dist_tss[dist_tss < 0]).argmin()
			gene_tss_down = gene[c_index[(dist_tss[dist_tss < 0][down_tss] == dist_tss).nonzero()[0]]][0]
			table.write('%s\t%s\t%s\t%s\t%s\t%s\t\t%d\t%d\t%f\t%f\n' % (gene_tss_down[3],gene_tss_down[2],gene_tss_down[4],gene_tss_down[5],gene_tss_down[6],gene_tss_down[7],peak.left[k],peak.right[k],peak.enrich1[k],peak.enrich2[k]))
		except ValueError:
			pass
		try:
			up_tes = dist_tes[dist_tes >= 0].argmin()
			gene_tes_up = gene[c_index[(dist_tes[dist_tes >= 0][up_tes] == dist_tes).nonzero()[0]]][0]
			table.write('%s\t%s\t%s\t%s\t%s\t%s\t\t%d\t%d\t%f\t%f\n' % (gene_tes_up[3],gene_tes_up[2],gene_tes_up[4],gene_tes_up[5],gene_tes_up[6],gene_tes_up[7],peak.left[k],peak.right[k],peak.enrich1[k],peak.enrich2[k]))
		except ValueError:
			pass
		try:
			down_tes = abs(dist_tes[dist_tes < 0]).argmin()
			gene_tes_down = gene[c_index[(dist_tes[dist_tes < 0][down_tes] == dist_tes).nonzero()[0]]][0]
			table.write('%s\t%s\t%s\t%s\t%s\t%s\t\t%d\t%d\t%f\t%f\n' % (gene_tes_down[3],gene_tes_down[2],gene_tes_down[4],gene_tes_down[5],gene_tes_down[6],gene_tes_down[7],peak.left[k],peak.right[k],peak.enrich1[k],peak.enrich2[k]))
		except ValueError:
			pass
			
	table.close()

def annotate_nooverlap(peak,ip,frac,gene,length = 2):
	"""Print a text file that gives a gene annotation for the overlapped peaks. It will print 
	the closest TSS and TES, up- and downstream of the middle of the peak window"""
	name = ip[0]
	for j in range(len(ip)-1):
		name = name + '_' + ip[j+1]
	
	print 'Printing the not overlapped gene annotation table'
	table = open('NoOverlap_%s_Top_%sPercent_Peaks.txt' % (name,str(frac)),'w')
	table.write('Gene Name\tGene ID\tChromosome\tStrand\tTSS\tTES\t\tPeak Left Edge\tPeak Right Edge\tReference IP Enrichment\n')
	
	# Create numpy arrays for faster computing
	chr = array(peak.chr)
	mid = array(peak.mid)
	
	for k in range(len(peak.chr)):
		# Only look at the genes on the same chromosome as the peak
		c = gene[:][:,4] == chr[k]
		
		# Output the index within the gene array of the appropriate genes
		c_index = c.nonzero()[0]
		
		# Return the TSS and TES values for the available genes
		tss = return_tss(gene[c,:])
		tes = return_tes(gene[c,:])
		
		# Calculate the distance to each TSS/TES from the peak midpoint
		dist_tss = tss - mid[k]
		dist_tes = tes - mid[k]
		
		# Find index of the closest upstream and downstream TSS and TES and return closest gene
		try:
			up_tss = dist_tss[dist_tss >= 0].argmin()
			gene_tss_up = gene[c_index[(dist_tss[dist_tss >= 0][up_tss] == dist_tss).nonzero()[0]]][0]
			table.write('%s\t%s\t%s\t%s\t%s\t%s\t\t%d\t%d\t%f\n' % (gene_tss_up[3],gene_tss_up[2],gene_tss_up[4],gene_tss_up[5],gene_tss_up[6],gene_tss_up[7],peak.left[k],peak.right[k],peak.enrich1[k]))
		except ValueError:
			pass
		try:
			down_tss = abs(dist_tss[dist_tss < 0]).argmin()
			gene_tss_down = gene[c_index[(dist_tss[dist_tss < 0][down_tss] == dist_tss).nonzero()[0]]][0]
			table.write('%s\t%s\t%s\t%s\t%s\t%s\t\t%d\t%d\t%f\n' % (gene_tss_down[3],gene_tss_down[2],gene_tss_down[4],gene_tss_down[5],gene_tss_down[6],gene_tss_down[7],peak.left[k],peak.right[k],peak.enrich1[k]))
		except ValueError:
			pass
		try:
			up_tes = dist_tes[dist_tes >= 0].argmin()
			gene_tes_up = gene[c_index[(dist_tes[dist_tes >= 0][up_tes] == dist_tes).nonzero()[0]]][0]
			table.write('%s\t%s\t%s\t%s\t%s\t%s\t\t%d\t%d\t%f\n' % (gene_tes_up[3],gene_tes_up[2],gene_tes_up[4],gene_tes_up[5],gene_tes_up[6],gene_tes_up[7],peak.left[k],peak.right[k],peak.enrich1[k]))
		except ValueError:
			pass
		try:
			down_tes = abs(dist_tes[dist_tes < 0]).argmin()
			gene_tes_down = gene[c_index[(dist_tes[dist_tes < 0][down_tes] == dist_tes).nonzero()[0]]][0]
			table.write('%s\t%s\t%s\t%s\t%s\t%s\t\t%d\t%d\t%f\n' % (gene_tes_down[3],gene_tes_down[2],gene_tes_down[4],gene_tes_down[5],gene_tes_down[6],gene_tes_down[7],peak.left[k],peak.right[k],peak.enrich1[k]))
		except ValueError:
			pass
			
	table.close()

def return_tss(g):
	"""Returns the TSS (in a numpy array) based off of what strand the gene finds itself"""
	tss = []
	for line in g:
		if line[5] == '+':
			tss.append(int(line[6]))
		else:
			tss.append(int(line[7]))
	return array(tss)
	
def return_tes(g):
	"""Returns the TES (in a numpy array) based off of what strand the gene finds itself"""
	tes = []
	for line in g:
		if line[5] == '+':
			tes.append(int(line[7]))
		else:
			tes.append(int(line[6]))
	return array(tes)