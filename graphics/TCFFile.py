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
		
		c = 29979245800.0 # speed of light in cm/s
		
		# set up the x axis - frequencies
		self.freq = [float(i)/N/dt/c for i in range(1,N/2)]
		
		# Here we'll do a little windowing to clean up the data
		self.fourier = []
		
		for data in self.data:
			data = map (lambda x,y: x * y, data, hamming(N))
			# now run an fft on the real data
			self.fourier.append(fft.rfft(data)[1:-1])

