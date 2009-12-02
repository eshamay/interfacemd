#!/usr/bin/python

import sys
import math
from ColumnDataFile import ColumnDataFile as CDF

import matplotlib.pyplot as plt
from pylab import *

def CheckNaN(x):
	if math.isnan(float(x)):
		return True
	else:
		return False

def SplitPairKey(key):
	key = key.replace('-',' ')
	key = key.replace('(',' ')
	key = key.replace(')',' ')
	key = key.split()
	return key

def PairEqual(left,right):
	if len(left) < 2 or len(right) < 2:
		return False
	elif left == right or (left[0] == right[1] and left[1] == right[0]):
		return True
	else:
		return False

def SortDictByKeys(dict):
	keys = dict.keys()
	keys.sort()
	return map(dict.get, keys)

def DataDictByPair(file,pair):
	ret = {}

	for key in file.keys():
		# check each key to see if it's the one we're looking for
		split_key = SplitPairKey(key)
		# Only work with the atom pair we're interested in
		if PairEqual(pair, split_key[:2]) and (float(split_key[2]) != 25.0 and float(split_key[2]) != 55.0):
			ret[float(split_key[2])] = file[key]

	return ret


class RDF2DPlotter:
	def __init__(self,file,pair,dims):	# The pair of atoms for the RDF (i.e. O O, or O H1, etc)
		self.file = CDF(file)
		self.pair = pair
		self.x = self.file['Distance']
		self.data = []

		self.InitData(pair)

		self.InitPlot(dims)
		if dims == 1:
			self.Plot1DData()
			self.ShowLegend()
		else:
			self.Plot2DData()

		self.PlotGraph()

	def InitPlot(self,dims):
		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		# some of the prelim stuff for our figure
		self.ax = self.fig.add_subplot(1,1,1)
		#self.ax.set_title(TITLE, size='x-large')
		self.cmBase = cm.Oranges   # choose a colormap to use for the plot


		xlabel(r'Inter-Atomic Distance / $\AA$', fontsize=28)
		xticks(fontsize=20)
		yticks(fontsize=20)

		if dims == 1:
			title('Radial Distribution Function', fontsize=36)
			ylabel(r'g(r)$_{'+self.pair[0]+self.pair[1]+'}$', fontsize=28)

		else:
			title(r'g(r)$_{'+self.pair[0]+self.pair[1]+'}$', fontsize=36)
			ylabel(r'Distance to Interface / $\AA$', fontsize=28)


	def Plot2DData(self):
		im = plt.imshow(self.data, cmap=self.cmBase, aspect='auto', extent=(0.0, 10.0, 30.0, 0.0))#, interpolation='bilinear')    

	def Plot1DData(self):
		for key in sorted(self.data_dict.iterkeys()):
			if 'nan' in self.data_dict[key]:
				print "BING!"
			self.ax.plot(self.x, self.data_dict[key], linewidth=4, label=str(key))

	def PlotGraph(self):

		draw()
		show()

	def ShowLegend(self):

		# set some legend properties.  All the code below is optional.  The
		# defaults are usually sensible but if you need more control, this
		# shows you how
		leg = plt.legend(loc='best', shadow=True)

		# the matplotlib.patches.Rectangle instance surrounding the legend
		frame = leg.get_frame()
		frame.set_facecolor('1.00')    # set the frame face color to light gray

		# matplotlib.text.Text instances
		for t in leg.get_texts():
			t.set_fontsize('x-large')    # the legend text fontsize

		# matplotlib.lines.Line2D instances
		for l in leg.get_lines():
			l.set_linewidth(4.0)  # the legend line width


	def InitData(self,pair):
		self.data_dict = DataDictByPair(self.file,pair)

		# before loading in the data, the keys need to be sorted by slab position
		for key in sorted(self.data_dict.iterkeys()):
			if reduce(lambda x,y: x and y, map(CheckNaN, self.data_dict[key])):
				del self.data_dict[key]
			else:
				self.data.append(self.data_dict[key])



rdf = RDF2DPlotter(sys.argv[1], sys.argv[2:4], int(sys.argv[4]))
