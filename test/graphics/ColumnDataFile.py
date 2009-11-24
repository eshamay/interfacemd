import numpy

class ColumnDataFile:

	def __init__(self,file):
		self.file = file

		self.ParseHeader()
		self.reader = numpy.loadtxt(self.file, skiprows=1)
		self.FormDataDict()

		# Output the file name of the data files we're working with
		print "Processing column data file: ", file
		# crunch through the density data

	def FormDataDict(self):
		self.data = {}
		for name in range(len(self.header)):
			self.data[self.header[name]] = self.reader[1:,name]

	def ParseHeader(self):
		self.reader = open(self.file)
		header = self.reader.next()
		header = header.split()
		self.header = header
		self.reader.close()
		return
	
	def __getitem__(self,item):
		return self.data[item]

	def __setitem__(self,k,v):
		self.data[k] = v
		return

	def fieldnames(self):
		return self.data.keys()

	def iteritems(self):
		return self.data.iteritems()
