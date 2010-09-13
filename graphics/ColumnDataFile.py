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
		if header != None:
			self.header = line1
		# Otherwise just use column indices as identifiers, and form a dict out of the data set
		else:
			self.header = range(len(line1))
			self.line1 = line1
			

 	def ParseData(self,file,header):
		self.data = [line.strip().split() for line in file.readlines()]
		if header == None:
			self.data.insert(0,self.line1)

		for i in self.data:
			if len(i) < 11:
				print i
		self.data = zip(*self.data)

		for col in range(len(self.data)):
			try:
				self.data[col] = [float(i) for i in self.data[col]]
			except ValueError:
				continue

		self.data = dict(zip(self.header, self.data))
	

	def __iter__(self):
		for data in self.data:
			yield data

	def __getitem__(self,index):
		return self.data[index]

	def keys(self):
		return self.data.keys()

	def __len__(self):
		return len(self.data)
