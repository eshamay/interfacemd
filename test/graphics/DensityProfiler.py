import csv
import numpy
import sys

from scipy import *
from scipy.optimize import leastsq
import scipy.io.array_import

import matplotlib.text
import matplotlib.pyplot as plt
import matplotlib.patches

from StatisticMachine import StatisticMachine as sm

class DensityProfiler:

	def __init__(self,files=[],avg=True):
		self.files = files

		# file i/o to get the data
		self.data = []
		self.fits = []

		self.avg = avg
		self.sm = sm()

		for file in self.files:
			print file
			# crunch through the density data
			d = self.DataDict(file)
			self.data.append (d)

			# perform analysis and find the fitting functions
			f = self.FitData(d)
			self.fits.append (f)

	def FitWater(self,data):
		x = data['position']
		h2o = data['h2o']

		# initial guess values
		p1 = 0.0
		p2 = 1.0
		p3 = 1.0
		p4 = 0.0
		x0 = 25.0
		d0 = 5.0
		x1 = 70.0
		d1 = 5.0

		if self.avg:
			p0 = array ([p1,p2,x0,d0])
			fit_func = self.sm.tanh_fit
		else:
			p0 = array ([p1,p2,p3,p4,x0,d0,x1,d1])
			fit_func = self.sm.double_tanh_fit

		residual=self.sm.FittingFunction(self.sm.residuals,function=fit_func)
		plsq = leastsq(residual, p0, args=(h2o,x), maxfev=100000)

		fit = fit_func(x,plsq[0])

		return (fit, p0)

	def FitIon(self,data,species):

		x = data['position']
		ion = list(data[species])

		# to set the correct position of the gaussian peak ('b')
		ion_max = max(ion)
		ion_max_index = ion.index(ion_max)
		b = x[ion_max_index]


		# initial guess values
		p1 = 0.0
		p2 = 1.0
		x0 = 0.0
		d0 = 5.0
		a = 1.0
		#b = -1.5
		c = 0.7

		# second interface
		x1 = 80.0
		d1 = 4.0
		e = 0.5
		f = 80.0
		g = 2.0

		if self.avg:
			p0 = array ([p1,p2,x0,d0,a,b,c])
			fit_func = self.sm.tanh_gaussian_fit
		else:
			p0 = array ([p1,p2,p2,p1,x0,d0,x1,d1,a,b,c,e,f,g])
			fit_func = self.sm.double_tanh_gaussian_fit

		residual = self.sm.FittingFunction(self.sm.residuals,function=fit_func)
		plsq = leastsq(residual, p0, args=(ion,x), maxfev=1000000)

		fit = fit_func(x,plsq[0])

		return (fit, p0)

	# Do fitting analysis
	def FitData(self,data):

		fit = {}
		# let's first get the fit data
		fit['position'] = data['position']
		fit['h2o'] = []
		fit['p_h2o'] = []
		fit['anion'] = []
		fit['p_anion'] = []
		fit['cation'] = []
		fit['p_cation'] = []

		# fit the water profile
		(fit['h2o'], fit['p_h2o']) = self.FitWater(data)

		if self.avg:
			# adjust the profile so the zero-point sits right on the zero-point of the fit
			# then we regather the position data
			for j in range(len(data['position'])):
				data['position'][j] = data['position'][j] - fit['p_h2o'][2]
				fit['position'][j] = data['position'][j]

			#recalculate the fit parameters with the new shifted zero-point and gather them, too
			(fit['h2o'], fit['p_h2o']) = self.FitWater(data)
			p0 = fit['p_h2o']

			# print out some data about the water fit
			print "Water Fit Data:"
			print "Z0 = % 8.3f\nd = % 8.3f\n90-10 = % 8.3f\npI = % 8.3f\n" % (p0[2], p0[3], p0[3]*2.197, p0[0])

		else:
			p0 = fit['p_h2o']
			print "Water fitting data\n\nZ0 = % 8.3f\nd = % 8.3f\n90-10 = % 8.3f\n\nZ1 = % 8.3f\nd = % 8.3f\n90-10 = % 8.3f\n" % (p0[4], p0[5], p0[5]*2.197, p0[6], p0[7], p0[7]*2.197)


		# Here we get the fit parameters for the cation and anion data (if we have those in the system)
		if self.avg and len(data['cation']) > 0:
			(fit['anion'], fit['p_anion']) = self.FitIon(data,'anion')
			(fit['cation'], fit['p_cation']) = self.FitIon(data,'cation')

			pA = fit['p_anion']
			pC = fit['p_cation']
			anion_max = max(fit['anion'])
			anion_max_index = list(fit['anion']).index(anion_max)
			anion_max_location = fit['position'][anion_max_index]
			fit['anion_max'] = anion_max_location
			cation_max = max(fit['cation'])
			cation_max_index = list(fit['cation']).index(cation_max)
			cation_max_location = fit['position'][cation_max_index]
			fit['cation_max'] = cation_max_location

			print "Anion\n\ta = % 8.3f\n\tb = % 8.3f\n\tc = % 8.3f\n\tPeak Center = % 8.3f" % (pA[4], pA[5], pA[6], anion_max_location)
			print "\n\tZ0 = % 8.3f\n\td = % 8.3f\n\t%s = % 8.3f\n" % (pA[2], pA[3], r'$\rho_I$', pA[0])
			print "Cation\n\ta = % 8.3f\n\tb = % 8.3f\n\tc = % 8.3f\n\tPeak Center = % 8.3f" % (pC[4], pC[5], pC[6], cation_max_location)
			print "\n\tZ0 = % 8.3f\n\td = % 8.3f\n\t%s = % 8.3f\n" % (pC[2], pC[3], r'$\rho_I$', pC[0])

			# lastly we'll use a statistics machine to calculate the goodness of the fits of these curves
			x_exp = data['position']
			x_fit = fit['position']

			w_r_square = self.sm.r_square((x_exp,data['h2o']), (x_fit,fit['h2o']))
			a_r_square = self.sm.r_square((x_exp,data['anion']), (x_fit,fit['anion']))
			c_r_square = self.sm.r_square((x_exp,data['cation']), (x_fit,fit['cation']))

			print "Goodness of fit\n\tWater = % 8.3f\n\tAnion = % 8.3f\n\tCation = % 8.3f\n\n" % (w_r_square, a_r_square, c_r_square)

		return fit

	def PlotData(self):

		# Set up the plot parameters (labels, size, limits, etc)
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		# the extents of the x-range
		xmin = 0.0
		xmax = 0.0

		for i in range(len(self.data)):
			# some of the prelim stuff for our figure
			ax = self.fig.add_subplot(len(self.data),1,i+1)

			ax.set_title(self.files[i], size='x-large')

			x = self.fits[i]['position']
			p0 = self.fits[i]['p_h2o']

			width_min = -p0[3]*2.197/2
			width_max = p0[3]*2.197/2

			width_extra = 5.0

			if self.avg == True:
				if xmin > width_min - width_extra:
					xmin = width_min - width_extra
				if xmax < width_max + width_extra:
					xmax = width_max + width_extra
			else:
				xmin = -5.0
				xmax = 100

			# first we'll plot the water data for each figure
			ax.plot(x, self.data[i]['h2o'], 'k:', linewidth=4, label=r'H$_2$O')
			# and also the water fit
			ax.plot(x, self.fits[i]['h2o'], 'k-', linewidth=3)

		# Then plot the anion and cation (if the system has any) and their fits (if we're doing an avg'd intface
			if len(self.data[i]['anion']) > 0:
				ax.plot(x, self.data[i]['anion'], 'r:', linewidth=4, label='Anion')
				ax.plot(x, self.data[i]['cation'], 'b:', linewidth=4, label='Cation')

			if self.avg:
				if len(self.data[i]['anion']) > 0:
					ax.plot(x, self.fits[i]['anion'], 'r-', linewidth=3)
					ax.plot(x, self.fits[i]['cation'], 'b-', linewidth=3)

					# some useful markers to show ion peak locations
					ax.axvline(self.fits[i]['anion_max'], color='r', linestyle='dotted', linewidth=2)
					ax.axvline(self.fits[i]['cation_max'], color='b', linestyle='dotted', linewidth=2)

				# on each plot shade the water region, and shade the interface region a different color
				### This is added to identify the water region with blue shading
				plt.axvspan(ax.get_xlim()[0],0.0, facecolor='b', alpha=0.2)
				### this should put a translucent rectangular area around the (90-10) interface region
				plt.axvspan(p0[2]+width_min,p0[2]+width_max, facecolor='g', alpha=0.2)

			ylim = 0.0
			if self.avg:
				for fit in self.fits:
					for key in fit.keys():
						if not key in ['anion', 'cation', 'h2o']:
							continue
						if key in ['anion', 'cation'] and len(fit[key]) == 0:
							continue
						l = fit[key]
						if max(l) > ylim:
							ylim = max(l) * 1.1
			else:
				ymax = max(self.data[i]['h2o'])
				ylim = ymax * 1.1

			ax.set_ylim([0.0,ylim])


			### this adds a box that clearly labels the width and 90-10 thickness of the interface from the fitting tanh
			'''
			gibbs = matplotlib.patches.Ellipse ((0, (p0[1]-p0[0])/2.0), 0.5, 0.1, alpha=0.0, fc="none")
			ax.annotate("Thickness = %5.3f\n\"90-10\" = %5.3f" % (p0[3], p0[3]*2.197), (0.0, (p0[1]+p0[0])/2.0),
				xytext=(100,40), textcoords='offset points', size=14,
				bbox=dict(boxstyle="round", fc=(1.0, 0.7, 0.7), ec=(1., .5, .5)),
				arrowprops=dict(arrowstyle="wedge,tail_width=1.",
					fc=(1.0, 0.7, 0.7), ec=(1., .5, .5),
					patchA=None,
					patchB=gibbs,
					relpos=(0.2, 0.8),
					connectionstyle="arc3,rad=-0.1"),
				)
			'''

			ax.set_axis_bgcolor('w')

			'''
			# set some legend properties.  All the code below is optional.  The
			# defaults are usually sensible but if you need more control, this
			# shows you how
			leg = plt.legend(loc='best', shadow=True)

			# the matplotlib.patches.Rectangle instance surrounding the legend
			frame = leg.get_frame()
			frame.set_facecolor('0.80')    # set the frame face color to light gray

			# matplotlib.text.Text instances
			for t in leg.get_texts():
				t.set_fontsize('medium')    # the legend text fontsize

			# matplotlib.lines.Line2D instances
			for l in leg.get_lines():
				l.set_linewidth(2.0)  # the legend line width
			'''

		#plt.savefig('density.avg.pdf', papertype='letter', dpi=600.0, edgecolor='w', orientation='landscape', format='pdf' )

		# set the xrange for each axis
		#for ax in self.fig.get_axes():
			ax.set_xlim([xmin,xmax])

		'''
		for fit in self.fits:
			print numpy.trapz(fit['h2o'],fit['position'],0.1)
			if len(fit['anion']) == 0:
				continue
			print numpy.trapz(fit['anion'],fit['position'],0.1)
			print numpy.trapz(fit['cation'],fit['position'],0.1)
		'''

		plt.show()


	# Should handle extracting all the data we want into a dictionary from the density profile data files
	def DataDict(self,filename):
		datareader = csv.reader(open(filename), dialect=csv.excel_tab)

		# the density data files (should) contain rows of the following:
		# 	position	h2o	ccl4	cation	anion
		# The position is just that - distance from the gibb's surface
		# The other data columns are densities... sort of. They're actually the densities divided by the molecular weights
		data = {}
		data['position'] = []
		data['h2o'] = []
		data['ccl4'] = []
		data['cation'] = []
		data['anion'] = []

		# process each row and extract all the relevant data for each species
		for row in datareader:
			row = row[0]
			row = row.strip()
			row = row.split()

			data['position'].append(float(row[0]))
			data['h2o'].append(float(row[1]))
			data['ccl4'].append(float(row[2]))
			if len(row) > 3:
				# ionic weights
				# na = 23.0
				# cl = 35.4527
				# no3 = 62.0
				# so4 = 96.06

				data['cation'].append(float(row[3]))
				data['anion'].append(float(row[4]))

		### this should make the fit closer by adjusting the data that is used in the plotting and calculations to only the section that we're interested in.
		index = 0
		if self.avg:
			index = data['position'].index(-12.0)
		else:
			index = data['position'].index(-5.0)

		data['position'] = data['position'][index:]
		data['h2o'] = data['h2o'][index:]
		data['ccl4'] = data['ccl4'][index:]
		if len(data['cation']) > 0:
			data['cation'] = data['cation'][index:]
			data['anion'] = data['anion'][index:]

		if self.avg:
			### Normalize all the data so that:
			#		1) bulk water is set to unity
			#		2) The area of each of the curves is equal

			# First we normalize all the areas to unity
			for key in ['h2o', 'anion', 'cation']:

				x = data['position']
				y = data[key]
				if key in ['anion', 'cation']:
					if len(data[key]) == 0:
						continue
				area = numpy.trapz(y,x,0.1)
				data[key] = data[key] / area

			# then we grab the water fit to find the bulk water density
			(fit, p0) = self.FitWater(data)
			for key in data.keys():
				if key == 'position':
					continue

				data[key] = data[key] / p0[0]

		return data
