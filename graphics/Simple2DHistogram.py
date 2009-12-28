#!/usr/bin/python
import sys
from scipy import *
import scipy.io.array_import
#from mpl_toolkits.axes_grid import AxesGrid


# plotting libs
import matplotlib.pyplot as plt 
from pylab import *

GIBBS = {'H2O':33.8755, 'NaCl':33.9399, 'NaNO3':31.5981, 'Na2SO4':40.50149}

class Histogram2DData:
  	def __init__ (self,file):
	  	self.name = file.split('/')[0]
		print self.name

		self.gibbs = GIBBS[self.name]
		print "gibbs = ", self.gibbs

		self.data = []
		datum = scipy.io.array_import.read_array(file)
  		self.maxPosition = self.gibbs+20.0
		print "maxPosition = ", self.maxPosition
  		self.size = int((self.maxPosition)/0.1)+200
		print "size = ", self.size

		for i in range(len(datum[0])):
		  	#self.data.append(datum[:,i])		# grab all data
		  	self.data.append(datum[:self.size,i])

	def Data(self):
	  	return self.data
	def Name(self):
	  	return self.name
	def Gibbs(self):
	  	return self.gibbs
	def MaxPosition(self):
		return self.maxPosition


class Histogram2D:
  	def __init__(self,files):
	  	self.data = []
		for file in files:
	  		self.data.append(Histogram2DData(file))

		self.InitPlot()

  		for i in range(len(self.data)):
		  	self.AddPlotData(i)

		self.PlotData()

	def InitPlot(self):
	  	self.fig = figure(num=1, facecolor='w', edgecolor='w', frameon=True)
		self.cm = cm.jet
		self.cax = axes([0.95, 0.33, 0.025, 0.33])

	def AddPlotData(self,num):
		self.fig.add_subplot(len(self.data)/2,2,num+1)

		gibbs = self.data[num].Gibbs()
  		max = self.data[num].MaxPosition()
  		print max

		ext = (-20.0-gibbs, max-gibbs, 1.0, 0.0)	# extents for 0-1 for the Normal vector
		if (num+1)%2:
			ext = (-20.0-gibbs, max-gibbs, 1.0, -1.0)	# extents for the Bisector vector

		im = plt.imshow(self.data[num].Data(), cmap=self.cm, aspect='auto', extent=ext, interpolation='bilinear')    

		# change the ticks on the colorbar (get rid of them?)
		cb = plt.colorbar(im, cax=self.cax) # grab the Colorbar instance
		for t in cb.ax.get_yticklabels():
  			t.set_fontsize(20)
  		cb.ax.set_yticklabels(())	


		xlim(-10.0,18.0)

	def PlotData(self):

	  			# set the xlimit to the interesting region
		# Pass go - collect $200
	  	draw()
		show()

histo = Histogram2D (sys.argv[1:])
