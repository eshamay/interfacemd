#!/usr/bin/python

import csv
import numpy
import sys
import os

from scipy import *

import matplotlib.text
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

class MoritaSFG:

	def __init__(self, files):

		# file i/o to get the data
		self.data = []
		for f in files:
			self.data.append(self.DataDict(f))
		
		#self.ScaleData()


	def ScaleData(self):
		for d in range(len(self.data)):
			chi_max = max(self.data[d]['chi'])
			scale = lambda a: a/chi_max
			self.data[d]['chi'] = map(scale, self.data[d]['chi'])

	def DataDict(self,file):
		d = {'x':[],'real':[],'imag':[],'comp':[],'chi':[]}

		data = loadtxt(file)
		d['file'] = file
		d['x'] = data[:,0]
		d['real'] = data[:,1]
		d['imag'] = data[:,2]

		# complex number representation
		for j in range(len(d['real'])):
			c = complex(d['real'][j],d['imag'][j])
			d['comp'].append(c)

		# the total chi-squared lineshapes
		for c in d['comp']:
			#chi.append(abs(c)*abs(c))
			d['chi'].append(real(c*c.conjugate()))

		return d

	def SetLegend(self):
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

	def PlotData(self):

		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
		colors = ['r','b','g','c','m','y','k']
		color = 0

		ax = self.fig.add_subplot(1,1,1)
		
		# plot all the plots
		for d in self.data:
			ax.plot(d['x'], d['chi'], colors[color]+'-', linewidth=3, label=d['file'])
			color = color + 1
			ax.set_xlabel('Frequency')
			ax.set_yticklabels([])

			labels = ax.get_xticklabels() + ax.get_yticklabels()
			for label in labels:
				label.set_size('x-large')

		self.SetLegend()
		plt.show()

d = MoritaSFG(sys.argv[1:])
d.PlotData()
