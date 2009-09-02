#!/usr/bin/python

import sys

import matplotlib.pyplot as plt
from numpy import *

def PlotData():
	
	# Set up the plot parameters (labels, size, limits, etc)
	fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

	# the extents of the x-range
	xmin = 2000.0
	xmax = 4000.0

	chi_ax = fig.add_subplot(3,1,1)
	re_ax = fig.add_subplot(3,1,2)
	im_ax = fig.add_subplot(3,1,3)

	chi_ax.plot(x, chi, 'k-', linewidth=2)
	re_ax.plot(x, re, 'k-', linewidth=2)
	im_ax.plot(x, im, 'k-', linewidth=2)

	plt.show()

	#C = (1-math.exp(-beta*hbar*freq))/3.0/hbar/c

	#print freq, C*freq*abs(fourier[i])**2
	#print freq, real(fourier[i])



## first open up the bond-length trajectory file
myfile = open(sys.argv[1])
#
data = []
# # read the data in and make sure it's converted into a nice (floating point) format for manipulation
for line in myfile.readlines():
	line = line.strip()
	#line = line.split()
	if line != "":
		data.append(float(line))

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
temp = 300.0
k = 1.3806504e-23
beta = 1.0/temp/k
hbar = 1.05457148e-34

scaling_factor = 0.9981	# a vibrational scaling factor for frequencies taken from http://srdata.nist.gov/cccbdb/vsf.asp

# set up the x axis - frequencies
x = []
for i in range(1,N/2):
	x.append(float(i)/N/dt/c)

# Here we'll do a little windowing to clean up the data
window = hamming(N)

for i in range(N):
	data[i] = data[i] * float(window[i])

# now run an fft on the real data
fourier = fft.rfft(data)

chi = []
re = []
im = []

for i in fourier[1:-1]:
	h = 1.0j*i
	chi.append((h*h.conjugate()).real)
	re.append(i.real)
	im.append(i.imag)

print len(im)
print len(re)
print len(chi)
print len(x)

PlotData()
