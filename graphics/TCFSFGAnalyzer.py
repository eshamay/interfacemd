import numpy
	
# speed of light in cm/s
c = 29979245800.0 
k = 1.3806504e-23
hbar = 1.05457148e-34
# and print the output of the fft magnitudes at each frequency. The conversion factor used turns the indices into freq's in cm^-1 from freqs in Hz (*10^15)
scaling_factor = 0.9981	# a vibrational scaling factor for frequencies taken from http://srdata.nist.gov/cccbdb/vsf.asp

# a class that can suck in a data set (single-column array of floats) and do some math manipulations on it
class TCFSFGAnalyzer:

  	def __init__(self,data):
	  	self.data = data	# incoming data should be a single array of floats
		self.numDataPoints = len(data)

	def CalcFFT(self,dt):	# dt = timestep length in seconds
		# x-axis for freq-domain data (in cm-1)
		self.freq = [float(i)/self.numDataPoints/dt/c for i in range(1,self.numDataPoints/2)]	

		# a little windowing (hamming) to clean up the data
		windowed = map (lambda x,y: x * y, self.data, numpy.hamming(self.numDataPoints))
		
		# now run an fft on the real data
		self.fft = numpy.fft.rfft(windowed)[1:-1]
		max_fft = max(self.fft)
		self.fft = [i/max_fft for i in self.fft]
		self.fft_re = [i.real for i in self.fft]
		self.fft_imag = [i.imag for i in self.fft]


	def CalcSFGChi(self,dt,temp):
	  	
		beta = 1.0/temp/k
		self.CalcFFT(dt)

		#self.C = [(1-numpy.exp(-beta*hbar*f))/3.0/hbar/c for f in self.freq]
		#chi = map(lambda const, f, val: const*f*abs(val)**2, self.C, self.freq, self.fft)
		chi = map(lambda x: x*1j, self.fft)

  		# normalization to unity
		max_chi = max(chi)
		self.chi = [i/max_chi for i in chi]

		return



	def Data(self):				# the original untouched TCF data
		return self.data
	def Freq(self):
	  	return self.freq		# the frequency axis
	def FFTNorm(self):
	  	return [abs(x) for x in self.fft]# the magnitude of the pure FFT of the TCF data
	def FFTNormSquared(self):
		return [(x*x.conjugate()).real for x in self.fft]
	def FFTComplex(self):
		return self.fft			# complex data from the FFT
	def FFTReal(self):
	  	return self.fft_re		# just the real portion
	def FFTImag(self):
	  	return self.fft_imag	# the imaginary fft data
	def ChiSquared(self):
		return self.chi			# the Chi(2)^2 from SFG data

# vim: tabstop=4 shiftwidth=4 softtabstop=4
