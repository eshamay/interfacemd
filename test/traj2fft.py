#!/usr/bin/python

import sys

import matplotlib.pyplot as plt
from numpy import *

def PlotData():
	
	# Set up the plot parameters (labels, size, limits, etc)
	fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

	chi_ax = fig.add_subplot(2,1,1)
	parts_ax = fig.add_subplot(2,1,2)

	chi_ax.plot(freq, chi, 'k-', linewidth=2)
	parts_ax.plot(freq, re, 'r-', linewidth=2)
	parts_ax.plot(freq, im, 'b-', linewidth=2)

	for ax in fig.get_axes():
	  ax.set_xlim(700.0, 4000.0)
	  #ax.set_ylim(0.0, 400.0)

	chi_ax.set_ylim(0.0,1.1)
	parts_ax.set_ylim(-800.0,800.0)
	plt.show()


## first open up the bond-length trajectory file
myfile = open(sys.argv[1])

# # read the data in and make sure it's converted into a nice (floating point) format for manipulation
data = [float(line.strip().split()[0]) for line in myfile.readlines()]

myfile.close()

N = len(data)	# number of samples/data points
dt = 0.75e-15 # sample size in seconds

# Here we take the magnitude squared of the fourier transform (as per V. Buch 2007)
#magsquare = []
#for i in range(len(fourier)):
	#magsquare.append ((abs(fourier[i]))**2)
# this produces an array of complex values corresponding to the intensities of the spectrum. Each index corresponds to a        frequency

# and print the output of the fft magnitudes at each frequency. The conversion factor used turns the indices into freq's in cm^-1 from freqs in Hz (*10^15)
c = 29979245800.0 # speed of light in cm/s
temp = 273.0
k = 1.3806504e-23
beta = 1.0/temp/k
hbar = 1.05457148e-34

scaling_factor = 0.9981	# a vibrational scaling factor for frequencies taken from http://srdata.nist.gov/cccbdb/vsf.asp

# set up the x axis - frequencies
freq = [float(i)/N/dt/c for i in range(1,N/2)]

# Here we'll do a little windowing to clean up the data
data = map (lambda x,y: x * y, data, hamming(N))

# now run an fft on the real data
fourier = fft.rfft(data)[1:-1]

#C = (1-math.exp(-beta*hbar*freq))/3.0/hbar/c

#print freq, C*freq*abs(fourier[i])**2	# chi^2
#print freq, real(fourier[i])

C = [(1-math.exp(-beta*hbar*f))/3.0/hbar/c for f in freq]
chi = map(lambda c, f, val: c*f*abs(val)**2, C, freq, fourier)
max_chi = max(chi)
chi = [c/max_chi for c in chi]
#chi = [C*freq*abs(f)**2 for freq,f in x,fourier]
#chi = [sqrt(i*i.conjugate()).real for i in fourier]
re = [i.real for i in fourier]
im = [i.imag for i in fourier]

PlotData()
