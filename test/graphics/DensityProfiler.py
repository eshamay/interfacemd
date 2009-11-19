from ColumnDataFile import ColumnDataFile as CDF
from DensityFitter import DensityFitter as DF

import matplotlib.pyplot as plt
import matplotlib.text
import matplotlib.patches

TITLE = 'Water in NaNO3'
XRANGE = [0.0,120.0]
class DensityProfiler:

	def __init__(self,files=[],avg=False):
		self.files = files
		if len(files) == 0:
			print "No files given. Try again with data filenames provided"
			raise NameError('No files given')

		# file i/o to get the data
		self.data = []

		self.avg = avg

		# crunch through the density data
		for file in self.files:
			self.data.append (CDF(file))

	def PlotData(self,cols=[],fit=True):

		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		colors = ['b','g','r','c','m','y','k','g','r','c','m']

		# the extents of the x-range
		#xmin = min(self.data[0][0])
		#xmax = max(self.data[0][0])

		ymax = 0.0

		for i in range(len(self.data)):
			# some of the prelim stuff for our figure
			ax = self.fig.add_subplot(len(self.data),1,i+1)
			ax.set_title(TITLE, size='x-large')

			x = self.data[i][0]
			size = self.data[i].NUMCOLUMNS

			#plot each column data-set against the x variable
			for j in range(size)[1:]:
				if len(cols) == 0:
					pass
				elif j not in cols:
					continue

				datum = self.data[i][j]
				datum = map(lambda x: x*18.01, datum)

				# sets the maximum of the graph
				dat_max = max(datum)
				if dat_max > ymax:
					ymax = dat_max

				# plots the data
				if j == 1:
					ax.plot(x, datum, colors[j-1]+':', linewidth=4, label=r'H$_2$O')
				else:
					ax.plot(x, datum, colors[j-1]+':', linewidth=4)

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


				# now take care of plotting the fitting functions
				if fit:
					if j == 1:
						self.PlotWaterFit(ax, x, datum, colors[j-1])		
			ylim = ymax * 1.1
			ax.set_ylim([0.0,ylim])
			ax.set_xlim(XRANGE)
			ax.set_axis_bgcolor('w')

		ax.set_xlabel(r'Slab Position $\AA$', size='x-large')
		ax.set_ylabel(r'$\rho_{H_2O}$ ($\frac{mg}{mL}$)', size='x-large')
		self.SetLegend()
		plt.show()

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



	def PlotWaterFit (self, ax, x, datum, color):
		
		d = DF()
		(fit, params) = d.FitWater(x, datum)
		ax.plot(x, fit, color+'-', linewidth=2, label=r'H$_2$O Fit')

		print "Bulk Density = %f\n" % (max(fit))



		### This is added to identify the water region with shading
		plt.axvspan(params[4]-params[5],params[4]+params[5], facecolor='b', alpha=0.2)
		plt.axvspan(params[6]-params[7],params[6]+params[7], facecolor='b', alpha=0.2)

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


