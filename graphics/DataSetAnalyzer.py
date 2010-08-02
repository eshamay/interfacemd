import numpy

# a class that can suck in a data set (single-column array of floats) and do some math manipulations on it
class DataSet:
  	def __init__(self,data):
	  	self.data = data	# incoming data should be a single array of floats
		self.numDataPoints = len(data)
	
	def CalcFFT(self,dt):	# dt = timestep/resolution
        # speed of light in cm/s
		c = 29979245800.0 
		# x-axis for freq-domain data (in cm-1)
		self.freq = [float(i)/self.numDataPoints/dt/c for i in range(1,self.numDataPoints/2)]	

		# a little windowing (hamming) to clean up the data
		self.data = map (lambda x,y: x * y, self.data, numpy.hamming(self.numDataPoints))

		# now run an fft on the real data
		self.fft = numpy.fft.rfft(self.data)[1:-1]
		self.fft_re = [i.real for i in self.fft]
		self.fft_imag = [i.imag for i in self.fft]

		return self.fft

	def CalcSFGChi(self,dt,temp):
		k = 1.3806504e-23
		beta = 1.0/temp/k
		hbar = 1.05457148e-34
		scaling_factor = 0.9981	# a vibrational scaling factor for frequencies taken from http://srdata.nist.gov/cccbdb/vsf.asp
	  	
		self.CalcFFT(dt)

		self.C = [(1-math.exp(-beta*hbar*f))/3.0/hbar/c for f in self.freq]
		chi = map(lambda c, f, val: c*f*abs(val)**2, self.C, self.freq, self.fft)

  		# normalization to unity
		max_chi = max(chi)
		self.chi = [c/max_chi for c in chi]

	def Data(self):
		return self.data
	def Freq(self):
	  	return self.freq
	def FFTNormSquared(self):
	  	return [abs(x) for x in self.fft]
		return norm_sq
	def FFTComplex(self):
		return self.fft
	def FFTReal(self):
	  	return self.fft_re
	def FFTImag(self):
	  	return self.fft_imag

# vim: tabstop=4 shiftwidth=4 softtabstop=4
