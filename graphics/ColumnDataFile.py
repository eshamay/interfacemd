'''
A data file with columns of data (or meta-data)
'''

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False


class ColumnDataFile:
	# need to supply which columns contain float data
  	def __init__(self,file, float_data_cols=None, header=None):
	  	self.filename = file
		fp = open(file)

		# Parse the header (1st line of the file) if there is one.
		# Otherwise just use column indices as identifiers, and form a dict out of the data set
		line1 = fp.readline().strip().split()
		self.data = [line.strip().split() for line in fp.readlines()]
		if header:
			self.header = line1
		else:
			self.data.insert(0,line1)
			
		self.data_array = zip(*self.data)

		if not header:
			self.header = range(len(self.data_array))

		if (float_data_cols):
			for col in float_data_cols:
				self.data_array[col] = [float(i) for i in self.data_array[col]]
		else:
			try:
				self.data_array = [[float(i) for i in col] for col in self.data_array]
			except ValueError:
				print "The data file has data that is not a number. The column numbers must be supplied to specify data columns"
				return

		self.data = dict(zip(self.header, self.data_array))
		fp.close()


	def __getitem__(self,item):
		return self.data[item]
	def iteritems(self):
		return self.data.iteritems()
	def keys(self):
		return self.data.keys()

