#!/usr/bin/python
import sys
from ColumnDataFile import ColumnDataFile
from DipPolAnalyzer import DipPolAnalyzer

from PlotUtility import *

import matplotlib.pyplot as plt

class MoritaSFG2002:

	def __init__(self, file, dt=0.75e-15):
		# first open up the bond-length trajectory file
	  	self.datafile = ColumnDataFile(file)
		self.dipoles = apply (zip, [self.datafile[i] for i in range(3)])
		self.polarizabilities = apply(zip, [self.datafile[i] for i in range(3,12)])
		self.dpa = DipPolAnalyzer(self.dipoles,self.polarizabilities)

	def PlotData(self):

		# Set up the plot parameters (labels, size, limits, etc)
		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		# the extents of the x-range
		xmin = 2000.0
		xmax = 4000.0

		# plotting out IR spectra
		ir = TCFSFGAnalyzer (self.dpa.IR_TCF())
		ir.CalcSFGChi(0.75e-15,273.0)
		axs = fig.add_subplot(1,1,1)
		axs.plot(ir.Freq()[1000:], ir.ChiSquared()[1000:], label='xyz IR')


		'''
		pol = [0,1,2]	# polarization combos
		for row in range(3):
			for col in range(2):

				pol = [pol[0],pol[2],pol[1]]
				axs = fig.add_subplot(3,2,(2*row)+(col+1))

				#for plotting out SFG data
				#sfg = TCFSFGAnalyzer (self.dpa.SFG_TCF(pol[1],pol[1],pol[0]))
				#sfg.CalcSFGChi(0.75e-15,273.0)
				#axs.plot(sfg.Freq(), sfg.ChiSquared(), label=str((pol[1],pol[1],pol[0])))

				#axs.plot(sfg.Freq()[1000:], [(x*x.conjugate()).real for x in sfg.FFTNormSquared()][1000:], label=str((pol[1],pol[1],pol[0])))
				axs.set_xlim(1000,4000)
				ShowLegend(axs)
				# swap the last 2 elements

				
			pol = pol[1:]+pol[:1]	# rotate to the next P polarization
		'''


sfg = MoritaSFG2002(sys.argv[1])

sfg.PlotData()
plt.show()
