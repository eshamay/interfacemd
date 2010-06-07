#!/usr/bin/python
import csv
import math
import numpy

# plotting libs
from pylab import *

# a class to load the data from a 2-d rdf file
class RDFFile:

  def __init__(self,filename):

	# read in the file's contents
	self.filename = filename
	
	# grab all the data
	d = [data for data in self.import_text(self.filename, ' ')]

	self.pairnames = d[0][1:]
  	self.data = self.parse_data(d[2:])

  def Plot(self, data, ax):

  	rdf_axis = [r[0] for r in data]
	datum = [r[1] for r in data]

	ax.plot(rdf_axis, datum, color='blue', linewidth=3, label=self.pairnames[0]+'--'+self.pairnames[1])
	return


  def import_text(self,filename,separator):
	for line in csv.reader(open(filename), delimiter=separator, skipinitialspace=True):
	  if line:
	    yield line

  #def parse_pairnames(self,text):
	## get all the atom pair names into a parsed list
	#print text.split()[1:]
	#return text.split()[1:]
	#pairs = map (lambda t: t.split(), text)
	# remove empty names
	#pairs = filter (lambda t: len(t)-1, pairs)
	#return pairs[0]

  def parse_params(self,text):
	text = map(lambda x: float(x), text[1:4]+text[5:8])
	return text[0:3], text[3:6]

  def parse_data(self,data):
	# Data coming in will be in string format - it needs to be converted to float
	# also we convert any NaN to 0.0
	row_convert = lambda row: [float(x) if not numpy.isnan(float(x)) else 0.0 for x in row]
  	return [row_convert(r) for r in data]


