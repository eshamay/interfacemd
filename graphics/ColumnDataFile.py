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
  	def __init__(self,file,header=None):
	  	self.filename = file
		fp = open(file)

		self.ParseHeader(fp,header)
		self.ParseData(fp,header)

		fp.close()


	def ParseHeader(self,file,header):
		# Parse the header (1st line of the file) if there is one.
		line1 = file.readline().strip().split()
		if header:
			self.header = line1
		# Otherwise just use column indices as identifiers, and form a dict out of the data set
		else:
			self.header = range(len(line1))
			self.line1 = line1
			
  

  	def ParseData(self,file,header):
		self.data = [line.strip().split() for line in file.readlines()]
		if not header:
			self.data.insert(0,self.line1)
		self.data = zip(*self.data)

		for col in range(len(self.data)):
		  	try:
				self.data[col] = [float(i) for i in self.data[col]]
			except ValueError:
				continue

#print "The data file has data that is not a number. The column numbers must be supplied to specify data columns"
#				return

#
#		if (float_data_cols):
#			for col in float_data_cols:
#				self.data_array[col] = [float(i) for i in self.data_array[col]]
#		else:
		self.data = dict(zip(self.header, self.data))


	def __iter__(self):
		for data in self.data:
			yield data

	def __getitem__(self,index):
		return self.data[index]

	def keys(self):
		return self.data.keys()

