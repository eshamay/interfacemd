#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
import sys
from scipy import *
import scipy.io.array_import
from scipy.optimize import leastsq

import matplotlib.pyplot as plt
import matplotlib.patches

from DensityProfiler import DensityProfiler

class OrderParameters:

	def __init__(self, system):


		self.data = self.DataDict('../OrderParameters/'+system+'.orderparams.avg.dat')
		density = DensityProfiler(['../Densities/'+system+'.density.avg.dat'])
		self.dens = density.data[0]
		self.fits = density.fits[0]
		self.names = {
					  'x':0,
					  'phi':1, 'theta':2, 'psi':3,
					  'sinphi':4, 'sintheta':5, 'sinpsi':6,
					  'sin2phi':7, 'sin2theta':8, 'sin2psi':9,
					  'sinsqphi':10, 'sinsqtheta':11, 'sinsqpsi':12,
					  'cosphi':13, 'costheta':14, 'cospsi':15,
					  'cos2phi':16, 'cos2theta':17, 'cos2psi':18,
					  'cossqphi':19, 'cossqtheta':20, 'cossqpsi':21,
					  'sinthetacos2phi':22, 'sinthetacos2psi':23,
					  'sinthetasqcos2phi':24, 'sinthetasqcos2psi':25,
					  'n':26
		}

	def DataDict(self, filename):

		data = []

		file = scipy.io.array_import.read_array(filename)

		for i in range(len(file[0])):
			data.append(file[:,i])

		return data

	def PlotData(self):

		d = self.data
		# add a figure for each file
		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=False)

		# okay, setting up the figure and subplots
		plt.subplots_adjust(wspace=0.2)

		S1ax = fig.add_subplot(2,1,1)
		S2ax = fig.add_subplot(2,1,2)
		#dens = fig.add_subplot(3,1,3)

		# some labels for the plot
		S1ax.set_title ('Water Order Parameters', size=25)
		S2ax.set_xlabel(r'Distance to Interface / $\AA$', size=20)
		S1ax.set_ylabel(r'S$_1=\frac{1}{2}\left<3\cos^2\theta-1\right>$', size=20)
		S2ax.set_ylabel(r'S$_2=\frac{\left<\sin^2\theta\cos2\psi\right>}{\left<\sin^2\theta\right>}$', size=20)

		xmin = -10
		xmax = 10
		x = d[self.names['x']]
		n = d[self.names['n']]

		s1 = d[self.names['cossqtheta']]
		# calculates S1 order parameter
		for i in range(len(s1)):
			s1[i] = 0.5*(3.0*s1[i]/n[i]-1.0)
		# plot the S1 curve
		S1ax.plot(x,s1,'k', linewidth=3)
		S1ax.set_ylim(-0.2,0.05)

		s2 = d[self.names['sinthetasqcos2psi']]
		s2_den = d[self.names['sinsqtheta']]
		for i in range(len(s2)):
			s2[i] = -s2[i] / s2_den[i]
		# plot the S2 curve
		S2ax.plot(x,s2,'k', linewidth=3)

		S2ax.set_ylim(-0.2,0.21)


		# shows the interface and plots vertical reference lines where the ions peak
		for ax in fig.get_axes():
			ax.axvspan(-self.fits['p_h2o'][3]*2.197/2, self.fits['p_h2o'][3]*2.197/2, facecolor='g', alpha=0.2)
			ax.axvspan(ax.get_xlim()[0],0.0, facecolor='b', alpha=0.2)
			if len(self.dens['anion']) > 0:
				ax.axvline(self.fits['anion_max'], color='r', linestyle='dotted', linewidth=3)
				ax.axvline(self.fits['cation_max'], color='b', linestyle='dotted', linewidth=3)
			ax.set_xlim([xmin,xmax])



		# setup a new axis for plotting the fits for reference
		ref_ax = S2ax.twinx()
		# plot the water fit for reference
		ref_ax.plot(self.fits['position'],self.fits['h2o'],'k:', linewidth=3, label=r'H$_2$O')
		# plot the ion fit lines for reference
		#if len(self.dens['anion']) > 0:
		#	ref_ax.plot(self.fits['position'],self.fits['anion'],'r:', linewidth=3, label='Anion')
		#	ref_ax.plot(self.fits['position'],self.fits['cation'],'b:', linewidth=3, label='Cation')

		ref_ax.set_xlim([xmin,xmax])

		'''
		# Set some legend properties
		#print tuple(coordinations + fits)
		#leg = ax.legend(coordinations + fits, 'best', shadow=True)
		leg = plt.legend(loc='best')

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

		plt.show()
'''
		S1ax.set_xticklabels([])
		S2ax.set_xticklabels([])
		S2ax.set_ylabel(r'S$_2$', size=20)
		rho.set_ylabel(r'Density / $\frac{g}{ml}$', size=20)

		for label in rho.get_xticklabels() + rho.get_yticklabels() + S1ax.get_yticklabels() + S2ax.get_yticklabels():
			label.set_fontsize(20)


		S1ax.set_ylim([-0.25, 0.05])
		S2ax.set_ylim([-0.3, 0.3])

		# set up a few vertical lines to show where the interface is
		S1ax.vlines(0.0, S1ax.get_ybound()[0], S1ax.get_ybound()[1], color='k', linestyles='dotted')
		S2ax.vlines(0.0, S2ax.get_ybound()[0], S2ax.get_ybound()[1], color='k', linestyles='dotted')

		S1ax.set_ylim([-0.25, 0.05])
		S2ax.set_ylim([-0.3, 0.3])
		rho.set_ylim([-.05,2.0])

		S1ax.set_xlim([-10,10])
		S2ax.set_xlim([-10,10])
		rho.set_xlim([-10,10])

		max_S2, max_index = max((s2, i) for i, s2 in enumerate(S2))
		min_S1, min_index = min((s1, i) for i, s1 in enumerate(S1))

		S1ax.annotate ("Min\n(% 5.3f,% 5.3f)" % (x[min_index], min_S1), (x[min_index], min_S1),
			xytext=(40,35), textcoords='offset points', size=20,
			arrowprops=dict(arrowstyle="->",
				connectionstyle="angle,angleA=180,angleB=110,rad=20"
				)
			)

		S2ax.annotate ("Max\n(% 5.3f,% 5.3f)" % (x[max_index], max_S2), (x[max_index], max_S2),
			xytext=(-120,-85), textcoords='offset points', size=20,
			arrowprops=dict(arrowstyle="->",
				connectionstyle="angle,angleA=0,angleB=-60,rad=30"
				)
			)


		# set some legend properties.  All the code below is optional.  The
		# defaults are usually sensible but if you need more control, this
		# shows you how
		#
		leg = rho.legend(('H2O', 'H2O-fit', 'CCl4'), 'best', shadow=True)
		# the matplotlib.patches.Rectangle instance surrounding the legend
		frame  = leg.get_frame()
		frame.set_facecolor('0.80')    # set the frame face color to light gray

		# matplotlib.text.Text instances
		for t in leg.get_texts():
			t.set_fontsize('medium')    # the legend text fontsize

		# matplotlib.lines.Line2D instances
		for l in leg.get_lines():
			l.set_linewidth(1.5)  # the legend line width

		#plt.savefig('density.dat.png', dpi=300.0, format='png' )

'''
