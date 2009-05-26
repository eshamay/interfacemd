import csv
import numpy
import sys

from scipy import *

import matplotlib.text
import matplotlib.pyplot as plt

class MoritaSFG:

	def __init__(self,file):
		self.file = file

		# file i/o to get the data
		self.data = loadtxt(file)
		self.x = self.data[:,0]
		self.real = self.data[:,1]
		self.imag = self.data[:,2]

		# complex number representation
		self.comp = []
		for i in range(len(self.real)):
			c = complex(self.real[i],self.imag[i])
			self.comp.append(c)

		self.mag = []
		for c in self.comp:
			self.mag.append(abs(c)*abs(c))

		self.real_mag = []
		for r in self.real:
			self.real_mag.append(abs(r)*abs(r))

		self.imag_mag = []
		for im in self.imag:
			self.imag_mag.append(abs(im)*abs(im))

	def PlotData(self):

		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		ax1 = self.fig.add_subplot(2,1,1)
		ax1.set_title('SFG Spectrum', size='x-large')
		ax2 = self.fig.add_subplot(2,1,2)
		ax2.set_xlabel(r'Frequency / $cm^{-1}$', size=20)


		# first we'll plot the water data for each figure
		ax2.plot(self.x, self.real, 'g', linewidth=1, label=r'Real')
		ax2.plot(self.x, self.imag, 'r', linewidth=1, label=r'Imag')
		#ax2.plot(self.x, self.real_mag, 'g:', linewidth=1, label=r'real magnitude')
		#ax2.plot(self.x, self.imag_mag, 'r:', linewidth=1, label=r'Imag magnitude')
		ax1.plot(self.x, self.mag, 'b', linewidth=1, label=r'Complex')
		ax1.set_axis_bgcolor('w')

		# set some legend properties.  All the code below is optional.  The
		# defaults are usually sensible but if you need more control, this
		# shows you how
		leg = plt.legend(loc='best', shadow=True)

		# the matplotlib.patches.Rectangle instance surrounding the legend
		frame = leg.get_frame()
		frame.set_facecolor('0.80')    # set the frame face color to light gray

		# matplotlib.text.Text instances
		for t in leg.get_texts():
			t.set_fontsize('medium')    # the legend text fontsize

		# matplotlib.lines.Line2D instances
		for l in leg.get_lines():
			l.set_linewidth(2.0)  # the legend line width

		plt.show()
