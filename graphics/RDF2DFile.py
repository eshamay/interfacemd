#!/usr/bin/python
import csv
import math
import numpy

# plotting libs
from pylab import *

# a class to load the data from a 2-d rdf file
class RDF2DFile:

  def __init__(self,filename):

	# read in the file's contents
	self.filename = filename
	
	# grab all the data
	d = [data for data in self.import_text(self.filename, ' ')]

	self.pairnames = self.parse_pairnames(d[0][1:])
  	(self.position_params, self.rdf_params) = self.parse_params(d[1])
  	self.data = self.parse_data(d[2:])
  	

  def colors(self):
	colors = ['r', 'b', 'g', 'o']
	for c in colors:
		yield c

  def colorgradient(self,num,hexcolor):
	dc = 0xffffff/num
	print dc
	for i in range(num):
		yield hex(dc*i)

	
  def Plot1D(self, data, ax, hilite=False):

	col = self.colorgradient(len(data[0]), 0xEF7300)
  	solid_col = self.colors()
  	rdf_axis = [r[0] for r in data]
  	for i in range(5,len(data[0])):

		#convert the hex color value to a string
		#c = '#'+str(hex(0x110a00*i)).replace('0x','').zfill(6)
	  	c = '#'+str(col.next().replace('0x','').zfill(6))
		print c
	  	
		width=2

		pos = i*self.position_params[2]+self.position_params[0]
		if hilite and pos in hilite:
			c = solid_col.next()
			width=4
			
		datum = [r[i] for r in data]
		ax.plot(rdf_axis, datum, color=c, linewidth=width, label=str(pos))

  def Plot2D(self, data):
	# create a figure
	fig = figure(1)

	im = plt.imshow(self.data, cmap=cm.jet, aspect='auto', extent=(self.position_params[0], self.position_params[1], self.rdf_params[1], self.rdf_params[0]), interpolation='bilinear')    

	axis([62.0, 75.0, 2.0, 15.0])	

	grid(True)

	# A little magic to create a colorbar - legend for the plot(s)
	cax = axes([0.95, 0.33, 0.025, 0.33])
	cb = plt.colorbar(im, cax=cax) # grab the Colorbar instance
	for t in cb.ax.get_yticklabels():
		t.set_fontsize(20)
	cb.ax.set_yticklabels(())	

	# Pass go - collect $200
	draw()
	show()

  def getBin(self, val):
	return int((val-self.position_params[0])/self.position_params[2])
	
  def import_text(self,filename,separator):
	for line in csv.reader(open(filename), delimiter=separator, skipinitialspace=True):
	  if line:
	    yield line

  def parse_pairnames(self,text):

	# get all the atom pair names into a parsed list
	pairs = text[1:]
	return pairs

  def parse_params(self,text):
	text = [float(t.split(':')[1]) for t in text if not t.isalpha()]
	return (text[0:3], text[3:6])

  def parse_data(self,data):
	# Data coming in will be in string format - it needs to be converted to float
	row_convert = lambda element: float(element)
  	nan_convert = lambda element: element if not numpy.isnan(element) else 0.0
  	set_convert = lambda row: map(nan_convert, map(row_convert, row))
  	data = map(set_convert, data)
  	return data
