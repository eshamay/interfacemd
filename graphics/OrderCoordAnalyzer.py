#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
import csv
import numpy
import sys

from DensityProfiler import DensityProfiler

import matplotlib.pyplot as plt
#import matplotlib.patches

class OrderCoord:
	def __init__(self,filename):
		self.file = filename

		self.data = self.DataDict(self.file + ".ordercoords.avg.dat")

		self.density = DensityProfiler([self.file + ".density.avg.dat"])
		self.fit_params = self.density.fits[0]['p_h2o']

	# Should handle extracting all the data we want into a dictionary (keys are the different coordination types) from the order-param-coordination files
	def DataDict(self,filename):
		datareader = csv.reader(open(filename), dialect=csv.excel_tab)

		# here we grab the names of each coordination type
		names = datareader.next()
		data = {}
		#initialize each coordination list/dict
		for key in names:
			if key == 'position':
				data[key] = []
				continue

			data[key] = {'S1':[], 'S2':[]}


		# next grab each data row
		for row in datareader:
			row = row[0].strip()
			row = row.split()

			# first grab the position
			data['position'].append(float(row[0]))

			# then for each coordination grab the S1/S2 pair and make a neat list of them all
			count = 1
			while count < len(names)-1:
				s1 = row[count*2-1]
				s2 = row[count*2]
				if s1 == 'nan' or s2 == 'nan':
					pair = (row[count*2-1],row[count*2])
				else:
					pair = (float(row[count*2-1]),float(row[count*2]))

				data[names[count]]['S1'].append(pair[0])
				data[names[count]]['S2'].append(pair[1])
				count = count + 1

		# delete empty entries because they screw things up
		for name in data.keys():
			if name == 'position':
				continue
			if len(data[name]['S1']) == 0:
				del data[name]

		# now tally the total order parameters (all of them added together), and just the ones for the free-oh
		size = len(data['position'])
		data['total'] = {'S1':[0.0]*size, 'S2':[0.0]*size}
		data['freeoh'] = {'S1':[0.0]*size, 'S2':[0.0]*size}
		total_count = [0] * size
		freeoh_count = [0] * size
		for name in data.keys():
			if name == 'position':
				continue

			for i in range(size):
				s1 = data[name]['S1'][i]
				s2 = data[name]['S2'][i]
				if s1 == 'nan' or s2 == 'nan':
					continue

				data['total']['S1'][i] = data['total']['S1'][i] + s1
				data['total']['S2'][i] = data['total']['S2'][i] + s2
				total_count[i] = total_count[i] + 1

				if name in ['O','OO','OOO','H','OH','OOH','OOOH']:
					data['freeoh']['S1'][i] = data['freeoh']['S1'][i] + s1
					data['freeoh']['S2'][i] = data['freeoh']['S2'][i] + s2
					freeoh_count[i] = freeoh_count[i] + 1

		# Do some averaging
		for i in range(size):
			if not float(total_count[i]) == 0:
				data['total']['S1'][i] = float(data['total']['S1'][i]) / float(total_count[i])
				data['total']['S2'][i] = float(data['total']['S2'][i])/ float(total_count[i])
			if not float(freeoh_count[i]) == 0:
				data['freeoh']['S1'][i] = float(data['freeoh']['S1'][i]) / float(freeoh_count[i])
				data['freeoh']['S2'][i] = float(data['freeoh']['S2'][i]) / float(freeoh_count[i])

		return data

	def PlotData(self):
		# Set up the plot parameters (labels, size, axes limits, etc)
		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=False)
		plt.subplots_adjust(wspace=0.5)
		axs = []
		coordinations = ['OH', 'OHH', 'OOH', 'OOHH', 'OOOHH', 'OOHHH']
		for i in range(2):

			ax = fig.add_subplot(2,1,i+1)
			axs.append(ax)
			ax.grid(False)
			ax.set_axis_bgcolor('w')
			# here we create S1 and S2 plots for each coordination type
			for coord in coordinations:

				if coord == 'position':
					continue

				# grab the position and order parameter data
				x = []
				S = []

				for j in range(len(self.data['position'])):
					if self.data[coord]['S1'] == 'nan' or self.data[coord]['S2'] == 'nan':
						continue

					x.append((self.data['position'][j]))
					S.append(self.data[coord]['S'+str(i+1)][j])

				ax.plot(x, S, linewidth=2, label=coord)

			if i == 0:
				ax.set_xticklabels([])
				ax.set_ylabel(r'S$_1$', size=20)
				ax.set_ylim([-0.25,0.1])
				plt.legend()
			if i == 1:
				ax.set_xlabel(r'Distance to Interface / $\AA$', size=20)
				ax.set_ylabel(r'S$_2$', size=20)
				ax.set_ylim([-0.5,0.5])

			### This is added to identify the water region with blue shading
			plt.axvspan(ax.get_xlim()[0],0.0, facecolor='b', alpha=0.2)
			### this should put a need translucent rectangular area around the interface for clear visualization
			plt.axvspan(-self.fit_params[3]*2.197/2,self.fit_params[3]*2.197/2, facecolor='g', alpha=0.2)

		'''
		# plot out a thick freeOH line on both axes
		# plot out the total order parameters for reference
		axs[0].plot(self.data['position'], self.data['freeoh']['S1'], 'k:', linewidth=3, label='Free-OH')
		axs[0].plot(self.data['position'], self.data['total']['S1'], 'k', linewidth=3, label='Total')

		axs[1].plot(self.data['position'], self.data['freeoh']['S2'], 'k:', linewidth=3, label='Free-OH')
		axs[1].plot(self.data['position'], self.data['total']['S2'], 'k', linewidth=3, label='Total')
		'''


		# overlay the water and ion density profiles on the lower plot for convenience
		fit = self.density.fits[0]
		x_fit = fit['position']
		h2o = fit['h2o']
		anion = fit['anion']
		cation = fit['cation']
		#axs[1].plot(x_fit, h2o, 'k:', linewidth=2, label=r'H$_2$O')

		# let's stick in the anion and cation locations for reference
		for ax in fig.get_axes():
			if len(anion) > 0:
				ax.axvline(self.density.fits[0]['anion_max'], color='r', linestyle='dotted', linewidth=2)
				ax.axvline(self.density.fits[0]['cation_max'], color='b', linestyle='dotted', linewidth=2)
				#ax.plot(x_fit, anion, 'r:', linewidth=2, label='Anion')
				#ax.plot(x_fit, cation, 'b:', linewidth=2, label='Cation')

			xmin = -self.fit_params[3] - 5.0
			xmax = self.fit_params[3] + 3.0
			ax.set_xlim([xmin,xmax])


		# Set some legend properties
		#leg = ax.legend(coordinations + fits, 'best', shadow=True)
		leg = plt.legend(loc='lower right')

		# the matplotlib.patches.Rectangle instance surrounding the legend
		frame = leg.get_frame()
		frame.set_facecolor('0.80')    # set the frame face color to light gray

		# matplotlib.text.Text instances
		for t in leg.get_texts():
			t.set_fontsize('medium')    # the legend text fontsize

		# matplotlib.lines.Line2D instances
		for l in leg.get_lines():
			l.set_linewidth(2.0)  # the legend line width

		plt.show()

	# Calculates theta and psi (based on S1 and S2) for a given location in the interface
	def CalcAngles(self,coord,location):

		i = self.data['position'].index(location)
		print location
		s1 = self.data[coord]['S1'][i]
		theta = numpy.arccos(numpy.sqrt((s1*2.0+1.0)/3.0))

		s2 = self.data[coord]['S2'][i]
		psi = numpy.arccos(s2)/2.0

		return (theta, psi)

	def OutputInfo(self):

		max_s2_oh = max(self.data['OH']['S2'][:230])
		max_s2_oh_index = self.data['OH']['S2'].index(max_s2_oh)
		max_s2_oh_position = self.data['position'][max_s2_oh_index]
		max_s2_ooh = max(self.data['OOH']['S2'][:235])
		max_s2_ooh_index = self.data['OOH']['S2'].index(max_s2_ooh)
		max_s2_ooh_position = self.data['position'][max_s2_ooh_index]

		(theta, psi) = self.CalcAngles('OH', max_s2_oh_position)
		theta = theta * 180.0 / numpy.pi
		psi = psi * 180.0 / numpy.pi
		print "OH data:\n\tmax s2 = %5.3f @ %5.3f\n\ttheta = %5.3f\n\tpsi = %5.3f\n" % (max_s2_oh, max_s2_oh_position, theta, psi)
		(theta, psi) = self.CalcAngles('OOH', max_s2_ooh_position)
		theta = theta * 180.0 / numpy.pi
		psi = psi * 180.0 / numpy.pi
		print "OOH data:\n\tmax s2 = %5.3f @ %5.3f\n\ttheta = %5.3f\n\tpsi = %5.3f\n" % (max_s2_ooh, max_s2_ooh_position, theta, psi)
