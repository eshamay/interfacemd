#!/usr/bin/python
import sys
from ColumnDataFile import ColumnDataFile
from DipPolAnalyzer import DipPolAnalyzer
from TCFSFGAnalyzer import TCFSFGAnalyzer

from PlotUtility import *

import matplotlib.pyplot as plt

class MoritaSFG2002:

	def __init__(self, file, dt=0.75e-15, temp=300.0):
		# first open up the bond-length trajectory file
		self.datafile = ColumnDataFile(file)
		self.dipoles = apply (zip, [self.datafile[i] for i in range(3)])
		self.polarizabilities = apply(zip, [self.datafile[i] for i in range(3,12)])

		self.dpa = DipPolAnalyzer(self.dipoles,self.polarizabilities)
		self.dt = dt
		self.temp = temp

	def PlotData(self):

		# Set up the plot parameters (labels, size, limits, etc)
		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
		fig2 = plt.figure(num=2, facecolor='w', edgecolor='w', frameon=True)

		# the extents of the x-range
		xmin = 2000.0
		xmax = 4000.0

		# plotting out IR spectra
		ir = TCFSFGAnalyzer (self.dpa.IR_TCF())
		ir.CalcSFGChi(self.dt,self.temp)
		axs = fig2.add_subplot(1,1,1)
		axs.plot(ir.Freq(), ir.ChiSquared(), label='xyz IR')
		axs.set_xlim(1000,4000)
		ShowLegend(axs)

		'''
		sfg = TCFSFGAnalyzer (self.dpa.SFG_TCF(0,0,2))
		sfg.CalcSFGChi(self.dt,300.0)
		axs = fig.add_subplot(2,1,2)
		axs.plot(sfg.Freq()[1000:], sfg.ChiSquared()[1000:], label='xyz SFG')
		'''


		pol = [0,1,2]	# polarization combos
		for row in range(3):

			axs = fig.add_subplot(3,1,row+1)

			# get both of the symmetric SSP data sets together
			# i.e. s1,s1,p & s2,s2,p
			sfg_tcf = self.dpa.SFG_TCF(pol[0],pol[1],pol[2])
			sfg = TCFSFGAnalyzer (sfg_tcf)
			sfg.CalcSFGChi(self.dt,self.temp)

			# now average them

			axs.plot(sfg.Freq(), sfg.ChiSquared(), label="P = "+str(pol[2]))

			axs.set_xlim(1000,4000)
			ShowLegend(axs)

			pol = pol[1:]+pol[:1]	# rotate to the next P polarization
			sys.stdout.flush()


sfg = MoritaSFG2002(sys.argv[1], 1.0e-15, 298.0)

sfg.PlotData()
plt.show()
