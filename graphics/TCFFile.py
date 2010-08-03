from numpy import *

class TCFFile:
	def __init__(self,file,temperature):
		## first open up the time-correlation data file
	  	self.filename = file
  		self.temperature = temperature

  		datafile = open(file)
  		self.ParseFileData(datafile)
		datafile.close()

	def ParseFileData(self,datafile):
		# read the data in and make sure it's converted into a nice (floating point) format for manipulation
		self.data = [line.strip().split() for line in datafile.readlines()]
		self.data = [[float(i) for i in row] for row in self.data]
		self.data = zip(*self.data)

		N = len(self.data[0])	# number of samples/data points
		dt = 0.75e-15 # sample size in seconds
		
		# Here we take the magnitude squared of the fourier transform (as per V. Buch 2007)
		# and print the output of the fft magnitudes at each frequency. The conversion factor used turns the indices into freq's in cm^-1 from freqs in Hz (*10^15)
		c = 29979245800.0 # speed of light in cm/s
		temperature = 273.0
		k = 1.3806504e-23
		beta = 1.0/self.temperature/k
		hbar = 1.05457148e-34
		
		scaling_factor = 0.9981	# a vibrational scaling factor for frequencies taken from http://srdata.nist.gov/cccbdb/vsf.asp
		
		# set up the x axis - frequencies
		self.freq = [float(i)/N/dt/c for i in range(1,N/2)]
		
		# Here we'll do a little windowing to clean up the data
		self.fourier = []
		
		for data in self.data:
			data = map (lambda x,y: x * y, data, hamming(N))
			# now run an fft on the real data
			self.fourier.append(fft.rfft(data)[1:-1])

			#C = [(1-math.exp(-beta*hbar*f))/3.0/hbar/c for f in self.freq]
			#chi = map(lambda c, f, val: c*abs(val)**2, self.C, self.freq, fourier)
			#chi = map(lambda c, val: c*abs(val)
			#self.chi.append([c/max_chi for c in chi])
			#chi = [C*freq*abs(f)**2 for freq,f in x,fourier]
			#chi = [sqrt(i*i.conjugate()).real for i in fourier]
			#self.re.append([i.real for i in fourier])
			#self.imag.append([i.imag for i in fourier])

