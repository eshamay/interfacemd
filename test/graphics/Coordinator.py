import csv
import numpy

import matplotlib.pyplot as plt
import matplotlib.patches
from DensityProfiler import DensityProfiler

class Coordinator:

	def __init__(self,files=[],norm=False):
		self.files = files
		self.data = []
		self.density = []
		self.normalize = False

		for file in files:
			d = self.DataDict(file + ".coordination.avg.dat")
			print file
			self.data.append(d)
			dens = DensityProfiler([file + ".density.avg.dat"])
			self.density.append(dens)

	def DataDict(self,file):
		data = {}

		datareader = csv.reader(open(file), dialect=csv.excel_tab)

		header = datareader.next()
		for i in range(len(header)):
			if header[i].isspace():
				del header[i]
				break

		for name in range(len(header)):
			data[header[name]] = []

		for row in datareader:
			row = row[0]
			row = row.strip()
			row = row.split()

			data['position'].append(float(row[0]))

			name = 1
			while name < len(header):
				data[header[name]].append(float(row[name]))
				name = name + 1

		return data

	# return the names of the coord type with the largest areas (under their curves). This might be useful if trying to plot only the types with the most substantial densities
	def TopAreas(self,data,num):

		coords = data.keys()
		del coords[coords.index('position')]
		areas = [0.0] * len(coords)

		# now calculate the area of each coords
		x = data['position']
		for i in range(len(coords)):
			area = numpy.trapz(data[coords[i]], x, dx=0.1)
			areas[i] = area

		# and make a list of the densest coords
		area_sort = sorted(areas)
		area_sort.reverse()

		tops = []
		for i in range(num):
			index = areas.index(area_sort[i])
			tops.append(coords[index])

		return tops

	def PlotData(self):

		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
		for i in range(len(self.data)):
			ax = fig.add_subplot(len(self.data),1,i+1)
			d = self.data[i]
			dens = self.density[i]
			x = d['position']
			ax.set_autoscale_on(True)

			for coord in self.TopAreas(d,10):
				if self.normalize == True:
					d[coord] = d[coord] / numpy.trapz(d[coord],x,dx=0.1)	# normalize by area under the curve
				ax.plot(x,d[coord], linewidth=2, label=coord)

			x_range = dens.fits[0]['p_h2o'][3]*2.197/2
			ax.axvspan(-x_range, x_range, facecolor='g', alpha=0.2)
			ax.axvspan(ax.get_xlim()[0],0.0, facecolor='b', alpha=0.2)
			ax.set_xlim([-x_range-6,x_range+6])

		leg = plt.legend(loc='best')
		# set some legend properties.  All the code below is optional.  The
		# defaults are usually sensible but if you need more control, this
		# shows you how

		# the matplotlib.patches.Rectangle instance surrounding the legend
		frame = leg.get_frame()
		frame.set_facecolor('0.80')    # set the frame face color to light gray

		# matplotlib.text.Text instances
		for t in leg.get_texts():
			t.set_fontsize('x-large')    # the legend text fontsize

		# matplotlib.lines.Line2D instances
		for l in leg.get_lines():
			l.set_linewidth(1.5)  # the legend line width

		plt.show()
