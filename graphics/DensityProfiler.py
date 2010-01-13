from ColumnDataFile import ColumnDataFile as CDF
from DensityFitter import DensityFitter as DF

import matplotlib.pyplot as plt
import matplotlib.text
import matplotlib.patches


ATOMS = ['O','C']
#ATOMS = ['O','C','NA','S']
#ATOMS = ['O','C','NA','N']
#ATOMS = ['O','C','NA','Cl']
DATA_LENGTH = 750
#ATOMS = ['O','C']
TITLE = r'The Aqueous Na$_2$SO$_4$ Solution - CCl$_4$ Interface'
XRANGE = [-9.0,17.5]
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
		molecular_weight = {'O':18.01, 'C':153.82, 'NA':22.99, 'N':62.00, 'S':96.06, 'Cl':35.45, 'SI':28.0855}
		scale = {'O':1.0, 'C':1.0, 'NA':10.0, 'N':5.0, 'S':5.0, 'Cl':5.0, 'SI':1.0}

		conversion = 1.0/(30.0*30.0*BIN_WIDTH * 1.0e-24 * 6.02e23)
		for k,v in self.data.iteritems():
			if k in ATOMS:
				self.data[k] = map(lambda x: x * molecular_weight[k] * scale[k] * conversion, v)

		return

	def PlotData(self,fit=True):

		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		# some of the prelim stuff for our figure
		ax = self.fig.add_subplot(1,1,1)
		#ax.set_title(TITLE, size='x-large')

		x = self.data['Position']

		# take care of plotting the fitting functions for the water
		if fit:
			self.PlotWaterFit(ax, x, self.data['O'], 'k-')

		x = self.ShiftAxis(x,self.shift)
		# plot the desired atomic densities
		for atom in ATOMS:

			datum = self.data[atom]
			if atom != "O" and atom != "C":
				print "yeehaw"
				print x[datum.index(max(datum[:DATA_LENGTH]))]
				ax.axvline(x[datum.index(max(datum[:DATA_LENGTH]))], color=LINESTYLES[atom][0], linestyle=':', linewidth=5)
				
				

			# sets the maximum of the graph
			#dat_max = max(datum)
			#if dat_max > YRANGE[1]:
				#YRANGE[1] = dat_max

			# plots the data
			ax.plot(x, datum, color=LINESTYLES[atom][0], linestyle='--', linewidth=5, label=LINESTYLES[atom][1])

		# Do some labeling of the water data curve extrema (max/min, peaks/troughs, etc) for various reasons
		#self.LabelWaterExtrema(ax,x,self.data['O'])

		ylim = YRANGE[1] * 1.1
		ax.set_ylim(YRANGE)
		ax.set_xlim(XRANGE)
		ax.set_axis_bgcolor('w')

		ax.set_xlabel(r'Distance to Interface ($\AA$)', size=40)
		ax.set_ylabel(r'$\rho_{H_2O}$ ($\frac{mg}{mL}$)', size=35)

		for a in ax.get_xticklabels() + ax.get_yticklabels():
			a.set_fontsize(35)

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
		fit_x = x[:DATA_LENGTH]
		(self.fit, self.params) = d.FitLowerWater(fit_x, datum[:DATA_LENGTH])
		shift = self.params["gibbs"]
		self.shift = shift
		width = self.params["width"]

		print "System: ", TITLE
		print "Left-side gibb's surface is located at: ", shift
		# plot the fit line for the water
		fit_x = self.ShiftAxis(fit_x,shift)
		ax.plot(fit_x, self.fit, linetype, linewidth=4, label=r'H$_2$O Fit')

		### This is added to identify the water region with shading
		plt.axvspan(-width/2.0,width/2.0, facecolor='b', alpha=0.2)

		# some text to distinguish the water regions
		# this adds a box that clearly labels the width and 90-10 thickness of the interface from the fitting tanh
		gibbs_1 = matplotlib.patches.Ellipse ((shift, max(datum)/2.0), 0.5, 0.1, alpha=0.0, fc="none")

		ax.annotate("Water/Decane Interface\nThickness = %5.3f\n\"90-10\" = %5.3f" % (width, width*2.197), (shift, max(datum)/2.0),
			    xytext=(-200,80), textcoords='offset points', size=14,
			    bbox=dict(boxstyle="round", fc=(1.0, 0.7, 0.7), ec=(1., .5, .5)),
			    arrowprops=dict(arrowstyle="wedge,tail_width=1.",
					    fc=(1.0, 0.7, 0.7), ec=(1., .5, .5),
					    patchA=None,
					    patchB=gibbs_1,
					    relpos=(0.2, 0.8),
					    connectionstyle="arc3,rad=-0.1"),
			    )

		'''
		gibbs_2 = matplotlib.patches.Ellipse ((params[6], max(datum)/2.0), 0.5, 0.1, alpha=0.0, fc="none")

		ax.annotate("Water Interface #2\nThickness = %5.3f\n\"90-10\" = %5.3f" % (params[7], params[7]*2.197), (params[6], max(datum)/2.0),
			    xytext=(100,80), textcoords='offset points', size=14,
			    bbox=dict(boxstyle="round", fc=(1.0, 0.7, 0.7), ec=(1., .5, .5)),
			    arrowprops=dict(arrowstyle="wedge,tail_width=1.",
					    fc=(1.0, 0.7, 0.7), ec=(1., .5, .5),
					    patchA=None,
					    patchB=gibbs_1,
					    relpos=(0.2, 0.8),
					    connectionstyle="arc3,rad=-0.1"),
			    )
		'''
	'''

	def PlotPDSWaterFit (self, ax, x, datum, color):
		
		d = DF()
		(fit, params) = d.FitWater(x, datum)
		ax.plot(x, fit, color+'-', linewidth=2, label=r'H$_2$O Fit')

		### This is added to identify the water region with shading
		plt.axvspan(params[4]-params[5],params[4]+params[5], facecolor='b', alpha=0.2)
		plt.axvspan(params[6]-params[7],params[6]+params[7], facecolor='#FD8700', alpha=0.2)

		# some text to distinguish the water regions
		# this adds a box that clearly labels the width and 90-10 thickness of the interface from the fitting tanh
		gibbs_1 = matplotlib.patches.Ellipse ((params[4], max(datum)/2.0), 0.5, 0.1, alpha=0.0, fc="none")

		ax.annotate("PDS Interface\nThickness = %5.3f\n\"90-10\" = %5.3f" % (params[5], params[5]*2.197), (params[4], max(datum)/2.0),
			    xytext=(-200,80), textcoords='offset points', size=14,
			    bbox=dict(boxstyle="round", fc=(1.0, 0.7, 0.7), ec=(1., .5, .5)),
			    arrowprops=dict(arrowstyle="wedge,tail_width=1.",
					    fc=(1.0, 0.7, 0.7), ec=(1., .5, .5),
					    patchA=None,
					    patchB=gibbs_1,
					    relpos=(0.2, 0.8),
					    connectionstyle="arc3,rad=-0.1"),
			    )

		gibbs_2 = matplotlib.patches.Ellipse ((params[6], max(datum)/2.0), 0.5, 0.1, alpha=0.0, fc="none")

		ax.annotate("Air Interface\nThickness = %5.3f\n\"90-10\" = %5.3f" % (params[7], params[7]*2.197), (params[6], max(datum)/2.0),
			    xytext=(100,80), textcoords='offset points', size=14,
			    bbox=dict(boxstyle="round", fc=(1.0, 0.7, 0.7), ec=(1., .5, .5)),
			    arrowprops=dict(arrowstyle="wedge,tail_width=1.",
					    fc=(1.0, 0.7, 0.7), ec=(1., .5, .5),
					    patchA=None,
					    patchB=gibbs_1,
					    relpos=(0.2, 0.8),
					    connectionstyle="arc3,rad=-0.1"),
			    )

		'''
		
def average(values):
	"""Computes the arithmetic mean of a list of numbers.
	>>> print average([20, 30, 70])
	40.0
	"""
	return sum(values, 0.0) / len(values)


