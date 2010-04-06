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
		self.file1 = CDF(file[0])
		self.file2 = CDF(file[1])
		self.file3 = CDF(file[2])
		self.pair = pair
		self.x = self.file1['Position']

		self.InitPlot()
  		column = self.PairDataColumn(pair)
  		self.data1 = self.file1[column]
  		self.data2 = self.file2[column]
  		self.data3 = self.file3[column]
		self.ax.plot(self.x, self.data1, color='k', linewidth=4, label=r'CCl$_4$')
		self.ax.plot(self.x, self.data2, color='b', linewidth=4, label='Hydrocarbon')
		self.ax.plot(self.x, self.data3, color='r', linewidth=4, label='Fluorocarbon')

		plt.xlim(0.0, 12.5)
  		plt.ylim(0.0, 1.5)

  		self.ShowLegend()

		self.PlotGraph()

  	def PairDataColumn(self,pair):
		return '('+pair[0]+'-'+pair[1]+')'

	def InitPlot(self):
		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		# some of the prelim stuff for our figure
		self.ax = self.fig.add_subplot(1,1,1)

		title('Radial Distribution Function', fontsize=36)
		xlabel(r'Inter-Atomic Distance / $\AA$', fontsize=36)
		#ylabel(r'g(r)$_{'+self.pair[0]+' - '+self.pair[1]+'}$', fontsize=36)
		ylabel(r'g(r)', fontsize=36)
		xticks(fontsize=20)
		yticks(fontsize=20)

		for a in self.ax.get_xticklabels() + self.ax.get_yticklabels():
			a.set_fontsize(40)

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
			t.set_fontsize(36)    # the legend text fontsize

		# matplotlib.lines.Line2D instances
		for l in leg.get_lines():
			l.set_linewidth(4.0)  # the legend line width



rdf = RDF2DPlotter(sys.argv[1:4], sys.argv[4:6])
#rdf = RDF2DPlotter(sys.argv[1], sys.argv[2:4])
