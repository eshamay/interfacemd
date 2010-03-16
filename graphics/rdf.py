#!/usr/bin/python

import sys
import math
from ColumnDataFile import ColumnDataFile as CDF

import matplotlib.pyplot as plt
from pylab import *

def SplitPairKey(key):
	key = key.replace('-',' ')
	key = key.replace('(',' ')
	key = key.replace(')',' ')
	key = key.split()
	return key

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
	def __init__(self,file,pair):	# The pair of atoms for the RDF (i.e. O O, or O H1, etc)
		self.file = CDF(file)
		self.pair = pair
		self.x = self.file['Position']

		self.InitPlot()
  		column = self.PairDataColumn(pair)
  		self.data = self.file[column]
		self.ax.plot(self.x, self.data, linewidth=4, label=column)

		self.PlotGraph()

  	def PairDataColumn(self,pair):
		return '('+pair[0]+'-'+pair[1]+')'

	def InitPlot(self):
		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		# some of the prelim stuff for our figure
		self.ax = self.fig.add_subplot(1,1,1)

		title('Radial Distribution Function', fontsize=36)
		xlabel(r'Inter-Atomic Distance / $\AA$', fontsize=28)
		ylabel(r'g(r)$_{'+self.pair[0]+' - '+self.pair[1]+'}$', fontsize=28)
		xticks(fontsize=20)
		yticks(fontsize=20)

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



rdf = RDF2DPlotter(sys.argv[1], sys.argv[2:4])
