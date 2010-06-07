#!/usr/bin/python
import numpy
#from scipy import *
import csv
import sys

import matplotlib.pyplot as plt

XRANGE = [-15.0, 15.0]
INTERFACES = {'H2O':33.8755, 'NaCl':33.9399, 'NaNO3':31.5981, 'Na2SO4':40.50149}
LEGEND_NAMES = {'H2O':r'H$_2$O', 'NaCl':'NaCl', 'NaNO3':r'NaNO$_3$', 'Na2SO4':r'Na$_2$SO$_4$'}

class ChargeFile:
	def __init__(self,filename):
		self.filename = filename
		self.x, self.y = numpy.loadtxt(filename, skiprows=1, unpack=True)
		self.name = filename.split('/')[0]
		self.label = LEGEND_NAMES[self.name]
		self.x = map(lambda x: x - INTERFACES[self.name], self.x)

class ChargePlotter:
	def __init__(self,files):
		self.chargefiles = [ChargeFile(f) for f in files]

		# set up the figure
		self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
		ax = self.fig.add_subplot(1,1,1)
		# plot each data set
		for c in self.chargefiles:
			ax.plot(c.x, c.y, linewidth=3, label=c.label)

		ax.set_xlim(XRANGE)
		ax.set_ylabel(r'$\bar{E}(z)\propto\frac{q}{r_z^2}$ / a.u.', fontsize=36)
		ax.set_xlabel('Distance to Interface / $\AA$', fontsize=36)
		labels = ax.get_xticklabels() + ax.get_yticklabels()
		map(lambda x: x.set_size(36), labels)

		self.SetLegend()

		plt.show()
			
	def SetLegend(self):
		# set some legend properties.  All the code below is optional.  The
		# defaults are usually sensible but if you need more control, this
		# shows you how
		leg = plt.legend(loc='best', shadow=True)

		# the matplotlib.patches.Rectangle instance surrounding the legend
		frame = leg.get_frame()
		frame.set_facecolor('1.00')    # set the frame face color to light gray

		# matplotlib.text.Text instances
		for t in leg.get_texts():
			t.set_fontsize('x-large')    # the legend text fontsize

		# matplotlib.lines.Line2D instances
		for l in leg.get_lines():
			l.set_linewidth(4.0)  # the legend line width


c = ChargePlotter(sys.argv[1:])
