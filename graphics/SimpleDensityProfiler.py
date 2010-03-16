from ColumnDataFile import ColumnDataFile as CDF
from DensityFitter import DensityFitter as DF

import matplotlib.pyplot as plt
import matplotlib.text
import matplotlib.patches

import re

ATOMS = ['O', 'C18', 'C17']
# Area of the interface
XSIZE = 33.168
YSIZE = 33.512
XRANGE = [-10.0,10.0]
YRANGE = [0.0,2.0]
# line color/style and labels
LINESTYLES = {'O':['black',r'H$_2$O'], 'C':['blue',r'CCl$_4$'], 'NA':['#088618',r'Na'], 'N':['#a10f05',r'NO$_3$'], 'SI':['#6805a1',r'SO$_4$'], 'S':['#6805a1',r'SO$_4$'], 'Cl':['#6e6f00',r'Cl']}
BIN_WIDTH = 0.1

class DensityProfiler:

	def __init__(self,file):
		self.file = file

		# crunch through the density data
		self.data = CDF(file)

		# and scale it to be in units of g/mL instead of number density
		self.ScaleDataToDensity()

	# scale each data set to adjust the number density to density in terms of g/mL
	def ScaleDataToDensity(self):
		molecular_weight = {'O':18.01, 'C':12.00, 'NA':22.99, 'N':62.00, 'S':96.06, 'Cl':35.45, 'SI':28.0855}
		# A scaling factor to further alter the lineshape vertical scale
		scale = {'O':1.0, 'C':1.0, 'NA':10.0, 'N':5.0, 'S':5.0, 'Cl':5.0, 'SI':1.0}

		conversion = 1.0/(XSIZE*YSIZE*BIN_WIDTH * 0.602)
		for k,v in self.data.iteritems():
			if k in ATOMS:
				c = self.StripNameDigits(k)
				self.data[k] = map(lambda x: x * molecular_weight[c] * scale[c] * conversion, v)

		return

	# strips digits from atom names
	def StripNameDigits(self, name):
		return [item for item in re.split('[0-9]', name)][0]

	def PlotData(self):

		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		# some of the prelim stuff for our figure
		ax = self.fig.add_subplot(1,1,1)

		x = self.data['Position']

		# take care of plotting the fitting functions for the water and then shifting the axis to the GDS location
		self.PlotWaterFit(ax, x, self.data['O'], 'k-')
		x = self.ShiftAxis(x,self.shift)

		# plot the desired atomic densities
		for atom in ATOMS:
			datum = self.data[atom]

			# strip the name down to the atom without any other digits
			atom = self.StripNameDigits(atom)

			# plot the data
			ax.plot(x, datum, color=LINESTYLES[atom][0], linestyle='--', linewidth=5, label=LINESTYLES[atom][1])

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

