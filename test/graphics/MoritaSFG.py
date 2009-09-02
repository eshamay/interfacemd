import csv
import numpy
import sys
import os

from scipy import *

import matplotlib.text
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

class MoritaSFG:

	def __init__(self):
		self.files = []
		self.names = []
		#for file in os.listdir('../'):
			#if file.find('sfg.1-') > -1:
				#self.files.append('../'+file)
#
				#print file
				#file = file.strip("sfg.1-")
				#self.names.append(file)
				#print file
		#self.files.append('../sfg.0.8-0.6.dat')
		#self.files.append('../sfg.dsw-1.dat')
		#self.files.append('/raid1/Analysis/H2O+CTC/sfg.dsw-2.dat')
		self.files.append('/raid1/Analysis/H2O+CTC/sfg.alpha-switch.dat')
		#self.files.append('../sfg.set1-cos-term.dat')
		#self.files.append('../sfg.set2.dat')
		#self.files.append('/raid1/Analysis/NaCl+CTC/sfg.set2.dat')
		#self.files.append('/raid1/Analysis/NaNO3+CTC/sfg.set2.dat')
		#self.files.append('/raid1/Analysis/NaNO3+CTC/sfg.set2.67.dat')
		#self.files.append('/raid1/Analysis/Na2SO4+CTC/sfg.set2.dat')
		#self.files.append('/raid1/Analysis/Na2SO4+CTC/sfg.set2-1.dat')
		#self.files.append('/raid1/Analysis/Na2SO4+CTC/sfg.set2.55.dat')
		#self.files.append('/raid1/Analysis/Na2SO4+CTC/sfg.set2.60.dat')
		#self.files.append('/raid1/Analysis/Na2SO4+CTC/sfg.set2.65.dat')
		#self.files.append('/raid1/Analysis/Na2SO4+CTC/sfg.set2.70.dat')
		#self.files.append('/raid1/Analysis/Na2SO4+CTC/sfg.set1.70.dat')
		#self.files.append('/raid1/Analysis/Na2SO4+CTC/sfg.set2.70.dat')
		#self.files.append('/raid1/Analysis/NaCl+CTC/sfg.set2.0.8-0.6.70.dat')
		#self.files.append('/raid1/Analysis/NaNO3+CTC/sfg.set2.0.8-0.6.60.dat')
		#self.files.append('/raid1/Analysis/Na2SO4+CTC/sfg.set2.0.8-0.6.70.dat')
		self.files.append('/raid1/Analysis/NaCl+CTC/sfg.alpha-switch.dat')
		self.files.append('/raid1/Analysis/NaNO3+CTC/sfg.alpha-switch.dat')
		self.files.append('/raid1/Analysis/Na2SO4+CTC/sfg.alpha-switch.dat')
		#self.files.append('/raid1/Analysis/NaCl+CTC/sfg.dsw-1.dat')
		#self.files.append('/raid1/Analysis/NaNO3+CTC/sfg.dsw-1.dat')
		#self.files.append('/raid1/Analysis/Na2SO4+CTC/sfg.dsw-1.dat')

		#self.files.append('../sfg.DSW6.dat')
		#self.files.append('../sfg.DSW7-OHH.dat')

		self.names = self.files

		# file i/o to get the data
		self.data = []
		self.x = []

		self.real = []
		self.imag = []

		for i in range(len(self.files)):
			self.data.append(loadtxt(self.files[i]))
			self.x.append(self.data[i][:,0])
			self.real.append(self.data[i][:,1])
			self.imag.append(self.data[i][:,2])

		# complex number representation
		self.comp = []
		scale = 0.0
		for i in range(len(self.files)):
			self.comp.append([])

			if i == 1:
				scale = -20.0
			if i == 2:
				scale = -55.0
			if i == 3:
				scale = 35.0

			for j in range(len(self.real[i])):
				c = complex(self.real[i][j]+scale,self.imag[i][j])
				self.comp[i].append(c)

		# the total chi-squared lineshapes
		self.chi = []
		for i in range(len(self.files)):
			self.chi.append([])
			for c in self.comp[i]:
				#self.chi[i].append(abs(c)*abs(c))
				self.chi[i].append(real(c*c.conjugate()))

		# normalizing by one method or another
		for i in range(len(self.files)):
			# find the area of the curve:
			area = numpy.trapz(self.chi[i], self.x[i], dx=0.1)
			# find the max point in the curve
			max_peak = max(self.chi[i])
			for j in range(len(self.chi[i])):
				# normalize by area
				#self.chi[i][j] = self.chi[i][j] / area
				# normalize by free-oh peak
				self.chi[i][j] = self.chi[i][j] / max_peak

	def PlotData(self):

		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		'''
		extra_axes = True

		if extra_axes == True:
			ax1 = self.fig.add_subplot(3,1,1)
		else:
			ax1 = self.fig.add_subplot(1,1,1)

		ax1.set_title('SFG Spectrum', size='x-large')
		ax1.set_ylabel(r'$\left|\chi^{(2)}\right|^2$', size='xx-large')

		if extra_axes == True:
			ax2 = self.fig.add_subplot(3,1,2)
			ax2.set_ylabel(r'Re $\chi^{(2)}$', size='xx-large')
			ax3 = self.fig.add_subplot(3,1,3)
			ax3.set_ylabel(r'Im $\chi^{(2)}$', size='xx-large')
		'''

		'''
		for ax in self.fig.get_axes():
			#ax.set_yticklabels([])
			#ax.axhline(xmin=2800.0/4000.0, xmax=3800.0/4000.0, color='k', linestyle='-')
			ax.yaxis.set_major_locator(MaxNLocator(1))
			ax.yaxis.grid(True)
			ax.xaxis.grid(False)
		'''

		axs = []
		colors = ['r-','g-','b-']
		titles = [r'CCl$_4$-NaCl', r'CCl$_4$-NaNO$_3$', r'CCl$_4$-Na$_2$SO$_4$', r'Salt Comparison']
		for i in range(4):
			axs.append(self.fig.add_subplot(2,2,i+1))
			if i == 3:
				for j in range(3):
					axs[i].plot(self.x[j+1], self.chi[j+1], colors[j], linewidth=3, label=titles[j])
			else:
				# first the water plot
				axs[i].plot(self.x[0], self.chi[0], 'k-', linewidth=3)
				# now the salt plots
				axs[i].plot(self.x[i+1], self.chi[i+1], colors[i], linewidth=3)
			axs[i].set_yticklabels([])
			labels = axs[i].get_xticklabels() + axs[i].get_yticklabels()
			for label in labels:
				label.set_size('x-large')
			plt.title(titles[i], size='xx-large')

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

		plt.show()
