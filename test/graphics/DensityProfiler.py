from ColumnDataFile import ColumnDataFile as CDF
from DensityFitter import DensityFitter as DF

import matplotlib.pyplot as plt
import matplotlib.text
import matplotlib.patches

ATOMS = ['O','C']
TITLE = 'Water in NaNO3'
XRANGE = [0.0,120.0]
YRANGE = [0.0,0.0]
# line color/style and labels
LINESTYLES = {'O':['k:',r'H$_2$O'], 'C':['b:',r'CCl$_4$'], 'NA':['#C2A606:',r'Na$^+$'], 'N':['#A10F05:',r'NO$_3^-$'], 'S':['#6805A1:',r'SO$_4^{-2}'], 'CL':['#CBB808:',r'Cl$^-$']}

class DensityProfiler:

	def __init__(self,file):
		self.file = file
		# crunch through the density data
		self.data = CDF(file)
		# and scale it to be in units of g/mL instead of number density
		self.ScaleDataToDensity()

	# scale each data set to adjust the number density to density in terms of g/mL
	def ScaleDataToDensity(self):
		molecular_weight = {'O':18.01, 'C':153.82, 'NA':22.99, 'N':62.00, 'S':96.06, 'CL':35.45}

		for k,v in self.data.iteritems():
			if k in ATOMS:
				self.data[k] = map(lambda x: x * molecular_weight[k] / 270.996, v)

		return

	def PlotData(self,fit=True):

		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		# some of the prelim stuff for our figure
		ax = self.fig.add_subplot(1,1,1)
		ax.set_title(TITLE, size='x-large')

		x = self.data['Position']

		# plot the desired atomic densities
		for atom in ATOMS:

			datum = self.data[atom]

			# sets the maximum of the graph
			dat_max = max(datum)
			if dat_max > YRANGE[1]:
				YRANGE[1] = dat_max

			# plots the data
			ax.plot(x, datum, LINESTYLES[atom][0], linewidth=4, label=LINESTYLES[atom][1])

		# Do some labeling of the water data curve
		self.LabelWaterExtrema(ax,x,self.data['O'])

		# now take care of plotting the fitting functions
		if fit:
			self.PlotWaterFit(ax, x, self.data['O'], 'k-')

		ylim = YRANGE[1] * 1.1
		ax.set_ylim([YRANGE[0],ylim])
		ax.set_xlim(XRANGE)
		ax.set_axis_bgcolor('w')

		ax.set_xlabel(r'Slab Position $\AA$', size='x-large')
		ax.set_ylabel(r'$\rho_{H_2O}$ ($\frac{mg}{mL}$)', size='x-large')
		self.SetLegend()
		plt.show()

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
		(fit, params) = d.FitWater(x, datum)
		# plot the fit line for the water
		ax.plot(x, fit, linetype, linewidth=2, label=r'H$_2$O Fit')

		print "Bulk Density = %f\n" % (average([params[1],params[2]]))

		### This is added to identify the water region with shading
		plt.axvspan(params[4]-params[5],params[4]+params[5], facecolor='b', alpha=0.2)
		plt.axvspan(params[6]-params[7],params[6]+params[7], facecolor='b', alpha=0.2)

		'''
		# some text to distinguish the water regions
		# this adds a box that clearly labels the width and 90-10 thickness of the interface from the fitting tanh
		gibbs_1 = matplotlib.patches.Ellipse ((params[4], max(datum)/2.0), 0.5, 0.1, alpha=0.0, fc="none")

		ax.annotate("Water/Decane Interface\nThickness = %5.3f\n\"90-10\" = %5.3f" % (params[5], params[5]*2.197), (params[4], max(datum)/2.0),
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


