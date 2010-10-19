'''
A class for dealing with column-data-files where every pair of 2 columns represents a histogram pair.
The first column of the pair is the data point, and the 2nd column is the population
'''
from ColumnDataFile import ColumnDataFile as CDF
from Utility import *
import matplotlib.pyplot as plt

class ColumnPairFile (CDF):

 	def __init__(self,file,header=None):
		CDF.__init__(self,file,header)
		#self.filename = file
		#filedata = open(file)
		#self.ParseData(filedata)
		#filedata.close()

	def ParseData(self,filedata,header):
		# separate out the column pairs for easy access (self.data[0] = first column pair, etc.)
		# self.data[n][0] = data value for column-pair n, and self.data[n][1] = population for the values
		self.data = [[float(f) for f in row.strip().split()] for row in filedata.readlines() if len(row.split()) > 2]
		self.data = [group(d,2) for d in self.data]
		self.data = zip(*self.data)
		self.data = [zip(*d) for d in self.data]

	def PlotData(self,axs,columnpair=0,scale1=1.0,scale2=1.0):
		try:
			f1 = [x*scale1 for x in self.data[columnpair][0]]
			f2 = [x*scale2 for x in self.data[columnpair][1]]
			header = self.header[columnpair]
			self.handle = axs.plot(f1, f2, label=self.filename)
		except IndexError:
			if columnpair >= len(self.data):
				print "Only %d columns in the column-pair data file, but column %d was requested!" % (len(self.data), columnpair+1)

	def ColumnPair(self,columnpair):
		return self.data[columnpair]

	def __len__(self):
		return len(self.data)
