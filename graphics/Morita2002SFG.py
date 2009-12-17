#!/usr/bin/python

import sys

import matplotlib.text
import matplotlib.pyplot as plt
from numpy import *

class MoritaSFG2002:

	def __init__(self, file):
		## first open up the bond-length trajectory file
		myfile = open(file)

		self.data = []
		file_data = loadtxt(myfile)
		self.data.append(file_data[:,0])
		self.data.append(file_data[:,1])
		self.data.append(file_data[:,2])

		myfile.close()

		self.FFT()
		#self.PlotData()

	def FFT(self):
		self.fft = []
		self.freq = []
		self.real = []
		self.im = []
		self.chi = []

		dt = 0.75e-15 # sample size in seconds
		c = 29979245800.0 # speed of light in cm/s

		for dat in self.data:
			N = len(dat)

			# apply a hamming window filter
			window = hamming(N)
			for i in range(N):
				dat[i] = dat[i] * float(window[i])

			# fourier transform the data
			fourier = fft.rfft(dat)

			# save it for later
			self.fft.append(fourier)

			# grab the real, imag, and magnitudes of each spectrum

			chi = []
			re = []
			im = []
			for f in fourier[1:-1]:
				chi.append((f*f.conjugate()).real)
				re.append(f.real)
				im.append(f.imag)
			self.real.append(re)
			self.im.append(im)
			self.chi.append(chi)

			freq = []
			# setup the frequency scale for each spectrum
			for i in range(1,N/2):
				freq.append(float(i)/N/dt/c)
			self.freq.append(freq)


	def PlotData(self):

		# Set up the plot parameters (labels, size, limits, etc)
		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		# the extents of the x-range
		xmin = 2000.0
		xmax = 4000.0

		axs = []
		axs.append(fig.add_subplot(3,1,1))
		axs.append(fig.add_subplot(3,1,2))
		axs.append(fig.add_subplot(3,1,3))

		for i in range(len(axs)):
			axs[i].plot(self.freq[i], self.chi[i], 'k-', linewidth=2)
			axs[i].plot(self.freq[i], self.real[i], 'r-', linewidth=1)
			axs[i].plot(self.freq[i], self.im[i], 'b-', linewidth=1)

		plt.show()


sfg = MoritaSFG2002(sys.argv[1])
