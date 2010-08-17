'''
A class for dealing with column-data-files where every pair of 2 columns represents a histogram pair.
The first column of the pair is the data point, and the 2nd column is the population
'''
from ColumnDataFile import ColumnDataFile as CDF
from Utility import *

class ColumnPairFile (CDF):


class ColumnPairFile:
  	def __init__(self,file):
	  	self.filename = file
		
	  	filedata = open(file)
		self.ParseData(filedata)
		filedata.close()

	def ParseData(self,filedata):
		# separate out the column pairs for easy access (self.data[0] = first column pair, etc.)
	  	# self.data[n][0] = data value for column-pair n, and self.data[n][1] = population for the values
	  	self.data = [[float(f) for f in row.strip().split()] for row in filedata.readlines() if len(row.split()) > 2]
		self.data = [group(d,2) for d in self.data]
		self.data = zip(*self.data)
  		self.data = [zip(*d) for d in self.data]

	def PlotData(self,axs,columnpair=0):
	  	axs.plot(self.data[columnpair][0], self.data[columnpair][1], label=self.filename)

