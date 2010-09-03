import numpy
import csv

from DensityFitter import DensityFitter as DF

import matplotlib.pyplot as plt

XSIZE = 3.0
YSIZE = 3.0
BIN_WIDTH = 0.1

XRANGE = [-20.0, 110.0]
YRANGE = [0.0, 25.0]
class DensityProfiler:

	def __init__(self,file,atoms):
		self.file = file
		self.atoms = atoms

		# crunch through the density data
		d = [data for data in self.import_text(self.file, ' ')]

		# grab the atom names
		self.keys = d[0][1:]

		# parse all the data file
  		self.data = self.parse_data(d[1:])
  	
		# and scale it to be in units of g/mL instead of number density
		self.ScaleDataToDensity()

	def import_text(self,filename,separator):
		for line in csv.reader(open(filename), delimiter=separator, skipinitialspace=True):
			if line:
				yield line

  	def parse_data(self,data):
		# Data coming in will be in string format - it needs to be converted to float
		row_convert = lambda element: float(element)
  		nan_convert = lambda element: element if not numpy.isnan(element) else 0.0
  		set_convert = lambda row: map(nan_convert, map(row_convert, row))
  		data = map(set_convert, data)
  		return data

	# scale each data set to adjust the number density to density in terms of g/mL
	def ScaleDataToDensity(self):
		molecular_weight = {'O':18.01, 'OW':18.01, 'HW':1.00, 'C':12.00, 'NA':22.99, 'N':62.00, 'S':96.06, 'Cl':35.45, 'SI':28.0855}
		# A scaling factor to further alter the lineshape vertical scale
		scale = {'O':1.0, 'OW':1.0, 'HW':1.0, 'C':1.0, 'NA':10.0, 'N':5.0, 'S':5.0, 'Cl':5.0, 'SI':1.0}

		# conversion for the differential volume and number density to g/mL
		conversion = 1.0/(XSIZE*YSIZE*BIN_WIDTH * 602)
  		for key in self.keys:
			name = key.strip('0123456789')
			ind = self.keys.index(key)
  			if name in molecular_weight:
					self.data[ind] = map(lambda x: x * molecular_weight[name] * scale[name] * conversion, self.data[ind])

		return

	def PlotData(self):

		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		# some of the prelim stuff for our figure
		ax = self.fig.add_subplot(1,1,1)

		x = [d[0] for d in self.data]

		# take care of plotting the fitting functions for the water and then shifting the axis to the GDS location
		#self.PlotWaterFit(ax, x, self.data['O'], 'k-')
		#x = self.ShiftAxis(x,self.shift)

		# plot the desired atomic densities
		for atom in [key for key in self.keys if key in self.atoms]:
			print atom

			ind = self.keys.index(atom)+1
			datum = [d[ind] for d in self.data]

			# plot the data
			ax.plot(x, datum, linestyle='--', linewidth=3, label=atom)

		# Do some labeling of the water data curve extrema (max/min, peaks/troughs, etc) for various reasons
		#self.LabelWaterExtrema(ax,x,self.data['O'])

		ax.set_ylim(YRANGE)
		ax.set_xlim(XRANGE)
		ax.set_axis_bgcolor('w')

		ax.set_xlabel(r'Distance to Water Interface ($\AA$)', size=40)
		ax.set_ylabel(r'$\rho$ ($\frac{mg}{mL}$)', size=45)

		for a in ax.get_xticklabels() + ax.get_yticklabels():
			a.set_fontsize(40)

		#self.SetLegend()
		plt.show()

	def ShiftAxis(self,axis,shift):
		return map(lambda p: p-shift, axis)

	def LabelWaterExtrema(self,ax,x,datum):
			# Show the location of the maximum water density peak and
			# the trough behind it
			peak = max(datum)
			peak_index = datum.index(peak)
			ax.axhline(y=peak, color='b', linestyle=':')
			#here's the trough
			trough_region = datum[peak_index:peak_index+80]
			trough = min(trough_region)
			trough_index = datum.index(trough)

			ax.axhline(y=trough, color='b', linestyle=':')

			print "max density = %f\ntrough density = %f\n" % (peak, trough)
			print "Max Location = %f\nTrough Location = %f\n" % (x[peak_index], x[trough_index])


	def SetLegend (self):
		# set some legend properties.  All the code below is optional.  The
		# defaults are usually sensible but if you need more control, this
		# shows you how
		leg = plt.legend(loc='best', shadow=True)

		# the matplotlib.patches.Rectangle instance surrounding the legend
		frame = leg.get_frame()
		frame.set_facecolor('0.80')    # set the frame face color to light gray

		# matplotlib.text.Text instances
		for t in leg.get_texts():
			t.set_fontsize('x-large')    # the legend text fontsize

		# matplotlib.lines.Line2D instances
		for l in leg.get_lines():
			l.set_linewidth(2.0)  # the legend line width

	def PlotWaterFit (self, ax, x, datum, linetype='k-'):
		
		d = DF()

  		# If only looking at a specific region of the data (i.e. only one interface and not both)
  		DATACUT = 1000

		#Calculate the fitting parameters and the fitted lineshape
		(self.fit, self.params) = d.FitLowerWater(x[:DATACUT], datum[:DATACUT])

		shift = self.params["gibbs"]
		self.shift = shift
		print "Left-side gibb's surface is located at: ", shift

		width = self.params["width"]
		print "surface width is: ", width

		# plot the fit line for the water
		x = self.ShiftAxis(x,shift)
		ax.plot(x[:DATACUT], self.fit, linetype, linewidth=4, label=r'H$_2$O Fit')

		### This is added to identify the water interface region with shading
		plt.axvspan(-width/2.0,width/2.0, facecolor='b', alpha=0.2)

