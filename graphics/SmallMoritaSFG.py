#!/usr/bin/python

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
	  	d = {}

		data = loadtxt(file)
		d['file'] = file
		d['x'] = data[:,0]
		d['real'] = data[:,1]
		d['imag'] = data[:,2]

		# complex number representation
		d['comp'] = [complex(r,i) for r,i in zip(d['real'], d['imag'])]

		# the total chi-squared lineshapes
		d['chi'] = [real(c*c.conjugate()) for c in d['comp']]

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
		#colors = ['r','b','g','c','m','y','k','o']
		color = 0

		ax = self.fig.add_subplot(1,1,1)
		#ax2 = self.fig.add_subplot(2,1,2)
		
		# plot all the plots
		for d in self.data:
			#ax.plot(d['x'], d['chi'], colors[color]+'-', linewidth=3, label=d['file'])
			ax.plot(d['x'], d['chi'], linewidth=3, label=d['file'])
			#ax2.plot(d['x'], d['real'], 'b-', linewidth=3, label='Re')
			#ax2.plot(d['x'], d['imag'], 'r-', linewidth=3, label='Im')
			color = color + 1

		ax.set_ylabel(r'$|\chi^{(2)}|^2$', fontsize=28)
			#ax2.set_ylabel(r'$\chi^{(2)}$', fontsize=28)
		ax.set_yticklabels([])
			#ax2.set_yticklabels([])

		ax.set_xlabel('Frequency', fontsize=28)
		#labels = ax.get_xticklabels() + ax.get_yticklabels() + ax2.get_xticklabels() + ax2.get_yticklabels() 
		labels = ax.get_xticklabels() + ax.get_yticklabels()
		for label in labels:
			label.set_size('x-large')

		self.SetLegend()
		plt.show()

d = MoritaSFG(sys.argv[1:])
d.PlotData()
