import csv

class ColumnDataFile:

	def __init__(self,file):
		self.filepath = file
		if self.filepath == "":
			print "No filename given!"

		# Here's where the data will be held

		# Output the file name of the data files we're working with
		print "Processing file: ", file
		# crunch through the density data
		self.DataDict()

	# Should handle extracting all the data we want into a dictionary from the density profile data files
	def DataDict(self):

		self.data = []
		datareader = csv.reader(open(self.filepath), dialect=csv.excel_tab)
		for row in datareader:
			row = row[0].split()

			if len(self.data) == 0:
				for i in range(len(row)):
					self.data.append([])

			for i in range(len(row)):
				self.data[i].append(self.num(row[i]))

		self.NUMCOLUMNS = len(self.data)
		self.NUMROWS = len(self.data[0])

	def num (self,s):
		try:
			return int(s)
		except ValueError:
			return float(s)

	# For retrieving data from the column data file
	def __getitem__(self,key):

		try:
			return self.data[key]
		except IndexError:
			print "ColumnDataFile.__getitem__ :: not enough columns in the data file!"
			print "there are ", len(self.data), " columns in the file,"
			print "column ", key+1, " was requested"
			return "Error"
