import csv
import numpy
import sys
import os

from scipy import *

import matplotlib.text
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

class SFGData:
	def __init__(self,file):
		self.filename = file
		self.systemName = file.split('/')[-2]

		self.data = loadtxt(file)

		self.freq = self.data[:,0]		# Frequency for each data point
		self.real = list(self.data[:,1])		# The real part
		self.imag = list(self.data[:,2])		# The imaginary part

		# The complex representation of the data
		self.comp = map(lambda r,i: complex(r,i), self.real, self.imag)

		# chi is the magnitude-squared of the complex lineshape
		self.chi = map(lambda c: real(c*c.conjugate()), self.comp)
	
	def Filename(self):
		return self.filename
	def SystemName(self):
		return self.systemName
	def Freq(self):
		return self.freq
	def Real(self):
		return self.real
	def Imag(self):
		return self.imag
	def Complex(self):
		return self.comp
	def Chi(self):
		return self.chi

	# returns the chi lineshape scaled by some value
	def ScaleChi_LowFreqs(self,scale):
		return map(lambda x: x * scale, self.chi)

	def NormalizeChiToMaximum(self):
		# find the max point in the curve
		max_peak = max(self.chi)
		return self.ScaleChi_LowFreqs(1.0/max_peak)

	def NormalizeChiToArea(self):
		# find the area of the curve:
		area = numpy.trapz(self.chi, self.x, dx=0.1)
		return self.ScaleChi_LowFreqs(1.0/area)
	



class MoritaSFG:

	def __init__(self,files):
		self.files = files

		# file i/o to get the data
		self.sfgdata = []
		self.h2o = SFGData(files[0])

		for file in self.files[1:]:
			self.sfgdata.append(SFGData(file))

		self.InitPlot()

		#titles = [r'CCl$_4$-NaCl', r'CCl$_4$-NaNO$_3$', r'CCl$_4$-Na$_2$SO$_4$', r'Salt Comparison']
		for dat,ax in zip(self.sfgdata,self.axs):

			self.fig.sca(ax)

			self.PlotChi(dat,ax)
			self.SetTicksAndLabels(dat,ax)
			self.ShowLegend(ax)

			plt.xlim(2800,3800)
		plt.show()

	def InitPlot(self):

		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		self.axs = []
		for i in range(len(self.sfgdata)):
			self.axs.append(self.fig.add_subplot(2,2,i+1))
			plt.title(self.sfgdata[i].SystemName(), size='xx-large')

		return

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



	# Plot out all the real and imaginary parts
	def PlotRealImag(self,dat,ax):
		ax.plot(dat.Freq(), dat.Real(), 'b-', linewidth=2, label='Real')
		ax.plot(dat.Freq(), dat.Imag(), 'g-', linewidth=2, label='Imaginary')

	def ScaleChi_LowFreqs(self,dat,scale,cut_freq):
	  	cut_ind = list(dat.Freq()).index(cut_freq)
		re = dat.Real()

		# cut the imaginary part into two and scale the lower one
		imag = dat.Imag()
		low_imag = imag[:cut_ind]
		low_imag = map(lambda im: im*scale, low_imag)		# scale the complex part somehow
		high_imag = imag[cut_ind:]

		new_imag = low_imag + high_imag

		comp = map(lambda r,i: complex(r,i), re, new_imag)
		chi = map(lambda c: real(c*c.conjugate()), comp)

		return chi

	def ScaleChi(self,dat,scale):
		re = dat.Real()

		# cut the imaginary part into two and scale the lower one
		imag = dat.Imag()
		new_imag = map(lambda im: im*scale, imag)		# scale the complex part somehow

		comp = map(lambda r,i: complex(r,i), re, new_imag)
		chi = map(lambda c: real(c*c.conjugate()), comp)

		return chi


	def NormalizeListToMax(self,dat):
		max_peak = max(dat)
		return map(lambda x: x/max_peak, dat)

	def MakeFreqList(self,min,max,numbins):
	  	df = (max-min)/float(numbins)
		return [min + x*df for x in range(numbins)]

	def StretchFreq(self,freq,new_min,cut_freq):
	  	cut_freq_id = list(freq).index(cut_freq)

	  	low = list(freq[:cut_freq_id])
		high = list(freq[cut_freq_id:])

		new_low = self.MakeFreqList(new_min, cut_freq, len(low))

  		return new_low + high


	# plot the chi-squared lineshape
	def PlotChi(self,dat,ax):
		cut_freq = 0.0
		low_freq = 0.0
		line_color = 'k-'

	  	scale = 1.0
		if dat.SystemName() == 'NaCl':
			cut_freq = 3600.0
			low_freq = 1500.0
		  	scale = 1.4
			line_color = 'g-'
		elif dat.SystemName() == 'NaNO3':
			cut_freq = 3550.0
			low_freq = 2000.0
		  	scale = 1.5
			line_color = 'r-'
		elif dat.SystemName() == 'Na2SO4':
			cut_freq = 3630.0
			low_freq = 1000.0
		  	scale = 0.7
			line_color = 'purple'


  		# scale the chi
		h2o = self.h2o.Chi()
		h2o = self.NormalizeListToMax(h2o)

		chi = self.ScaleChi (dat, scale)
		chi = self.NormalizeListToMax(chi)

  		# scale the frequencies
  		dat_freq = self.StretchFreq (dat.Freq(), low_freq, cut_freq)
  		#dat_freq = dat.Freq()
  		h2o_freq = self.StretchFreq (self.h2o.Freq(), 2000, 3600.0)
  		#h2o_freq = self.h2o.Freq()

		#ax.plot(dat.Freq(), dat.Chi(), 'b:', linewidth=4, label=dat.SystemName()) #label=r'$|\chi|^{2}$')
		ax.plot(dat_freq, chi, line_color, linewidth=4, label=dat.SystemName()) #label=r'$|\chi|^{2}$')
		ax.plot(h2o_freq, h2o, 'k:', linewidth=4, label=r'H$_2$O')



	def SetTicksAndLabels(self,dat,ax):
		#axs[i].set_yticklabels([])
		labels = ax.get_xticklabels() + ax.get_yticklabels()
		for label in labels:
			label.set_size('x-large')

	def ShowLegend(self,ax):

		# set some legend properties.  All the code below is optional.  The
		# defaults are usually sensible but if you need more control, this
		# shows you how
		leg = plt.legend(loc='best', shadow=True, fancybox=True)

		# the matplotlib.patches.Rectangle instance surrounding the legend
		frame = leg.get_frame()
		frame.set_facecolor('1.00')    # set the frame face color to light gray

		# matplotlib.text.Text instances
		for t in leg.get_texts():
			t.set_fontsize('x-large')    # the legend text fontsize

		# matplotlib.lines.Line2D instances
		for l in leg.get_lines():
			l.set_linewidth(4.0)  # the legend line width

