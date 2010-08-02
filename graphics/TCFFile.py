from numpy import *

class TCFFile:
	def __init__(self,file,temp):
		## first open up the time-correlation data file
	  	self.filename = file
  		datafile = open(file)
  		self.temp = temp
  		self.ParseFileData(datafile)
		datafile.close()

	def ParseFileData(self,datafile):
		# read the data in and make sure it's converted into a nice (floating point) format for manipulation
		self.data = [float(line.strip().split()[0]) for line in datafile.readlines()]
		
		N = len(self.data)	# number of samples/data points
		dt = 0.75e-15 # sample size in seconds
		
		# Here we take the magnitude squared of the fourier transform (as per V. Buch 2007)
		#magsquare = []
		#for i in range(len(fourier)):
			#magsquare.append ((abs(fourier[i]))**2)
		# this produces an array of complex values corresponding to the intensities of the spectrum. Each index corresponds to a        frequency
		
		# and print the output of the fft magnitudes at each frequency. The conversion factor used turns the indices into freq's in cm^-1 from freqs in Hz (*10^15)
		c = 29979245800.0 # speed of light in cm/s
		#temp = 273.0
		k = 1.3806504e-23
		beta = 1.0/self.temp/k
		hbar = 1.05457148e-34
		
		scaling_factor = 0.9981	# a vibrational scaling factor for frequencies taken from http://srdata.nist.gov/cccbdb/vsf.asp
		
		# set up the x axis - frequencies
		self.freq = [float(i)/N/dt/c for i in range(1,N/2)]
		
		# Here we'll do a little windowing to clean up the data
		self.data = map (lambda x,y: x * y, self.data, hamming(N))
		
		# now run an fft on the real data
		fourier = fft.rfft(self.data)[1:-1]
		
		#C = (1-math.exp(-beta*hbar*freq))/3.0/hbar/c
		
		#print freq, C*freq*abs(fourier[i])**2	# chi^2
		#print freq, real(fourier[i])
		
		self.C = [(1-math.exp(-beta*hbar*f))/3.0/hbar/c for f in self.freq]
		chi = map(lambda c, f, val: c*f*abs(val)**2, self.C, self.freq, fourier)
		max_chi = max(chi)
		self.chi = [c/max_chi for c in chi]
		#chi = [C*freq*abs(f)**2 for freq,f in x,fourier]
		#chi = [sqrt(i*i.conjugate()).real for i in fourier]
		self.re = [i.real for i in fourier]
		self.im = [i.imag for i in fourier]

