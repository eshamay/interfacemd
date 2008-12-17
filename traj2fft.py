#!/usr/bin/python

import sys
from numpy import *

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

# Here we'll do a little windowing to clean up the data
window = hamming(N)

for i in range(N):
	data[i] = data[i] * float(window[i])

# now run an fft on the real data
fourier = fft.rfft(data)

# Here we take the magnitude squared of the fourier transform (as per V. Buch 2007)
magsquare = []
for i in range(len(fourier)):
	magsquare.append ((abs(fourier[i]))**2)
# this produces an array of complex values corresponding to the intensities of the spectrum. Each index corresponds to a        frequency

# and print the output of the fft magnitudes at each frequency. The conversion factor used turns the indices into freq's in cm^-1 from freqs in Hz (*10^15)
c = 29979245800.0 # speed of light in cm/s
temp = 300.0
k = 1.3806504e-23
beta = 1.0/temp/k
hbar = 1.05457148e-34

scaling_factor = 0.9981	# a vibrational scaling factor for frequencies taken from http://srdata.nist.gov/cccbdb/vsf.asp
for i in range(1,N/2):
	freq = float(i)/N/dt/c
	C = (1-math.exp(-beta*hbar*freq))/3.0/hbar/c

	print freq, C*freq*abs(fourier[i])**2
