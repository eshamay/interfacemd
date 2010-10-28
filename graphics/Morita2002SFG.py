#!/usr/bin/python

# This SFG calculator creates both IR and SFG spectra from dipole/polarizability files
# The SFG spectra are all SSP polarized, but that can be changed in the DipPol analyzer

# PrintData prints out 5 columns for frequency and spectral intensity:
#		frequency, IR, SFG_X, SFG_Y, SFG_Z

# PlotData creates 2 figures, one for IR and one for SFG. The SFG figure has 3 axes for X,Y, and Z polarization choices for dipole vector component

import sys
from ColumnDataFile import ColumnDataFile
from DipPolAnalyzer import DipPolAnalyzer
from TCFSFGAnalyzer import TCFSFGAnalyzer

#from PlotUtility import *

#import matplotlib.pyplot as plt
# the extents of the x-range
xmin = 1000.0
xmax = 15000.0


class MoritaSFG2002:

	def __init__(self, file, dt=0.75-15, temp=300.0):
		# first open up the bond-length trajectory file
		self.datafile = ColumnDataFile(file)
		self.dipoles = apply (zip, [self.datafile[i] for i in range(3)])
		self.polarizabilities = apply(zip, [self.datafile[i] for i in range(3,12)])

		self.dpa = DipPolAnalyzer(self.dipoles,self.polarizabilities)
		self.dt = dt
		self.temp = temp

	def PrintData(self):

		# plotting out IR spectra
		ir = TCFSFGAnalyzer (self.dpa.IR_TCF())
		ir.CalcSFGChi(self.dt,self.temp)
		ir_x = ir.Freq()
		ir_y = ir.FFTNormSquared()

		sfg_data = []

		pol = [0,1,2]	# polarization combos
		for row in range(3):

			# get both of the symmetric SSP data sets together
			# i.e. s1,s1,p & s2,s2,p
			sfg_tcf = self.dpa.SFG_TCF(pol[0],pol[1],pol[2])
			sfg = TCFSFGAnalyzer (sfg_tcf)
			sfg.CalcSFGChi(self.dt,self.temp)

			sfg_data.append(sfg.SFG())

			pol = pol[1:]+pol[:1]	# rotate to the next P polarization

		# limit the output to be within the frequency extents listed above (xmin, xmax)
		freq_min_index = ir_x.index([f for f in ir_x if f > xmin][0])
		freq_max_index = ir_x.index([f for f in ir_x if f < xmax][-1])
		# now get the data all set up
		data = zip(ir_x, ir_y, sfg_data[0], sfg_data[1], sfg_data[2])

		for d in data[freq_min_index:freq_max_index]:
			for i in d:
				print "%12.6e " % (i),
			print



'''
	def PlotData(self):

		# Set up the plot parameters (labels, size, limits, etc)
		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
		fig2 = plt.figure(num=2, facecolor='w', edgecolor='w', frameon=True)

		# the extents of the x-range
		xmin = 0.0
		xmax = 10000.0

		# plotting out IR spectra
		ir = TCFSFGAnalyzer (self.dpa.IR_TCF())
		ir.CalcSFGChi(self.dt,self.temp)
		axs = fig2.add_subplot(1,1,1)
		axs.plot(ir.Freq(), ir.ChiSquared(), label='xyz IR')
		axs.set_xlim(1000,4000)
		ShowLegend(axs)

		pol = [0,1,2]	# polarization combos
		for row in range(3):

			axs = fig.add_subplot(3,1,row+1)

			# get both of the symmetric SSP data sets together
			# i.e. s1,s1,p & s2,s2,p
			sfg_tcf = self.dpa.SFG_TCF(pol[0],pol[1],pol[2])
			sfg = TCFSFGAnalyzer (sfg_tcf)
			sfg.CalcSFGChi(self.dt,self.temp)

			# now average them

			axs.plot(sfg.Freq(), sfg.SFG(), label="P = "+str(pol[2]))

			axs.set_xlim(1000,4000)
			ShowLegend(axs)

			pol = pol[1:]+pol[:1]	# rotate to the next P polarization
			sys.stdout.flush()
'''


sfg = MoritaSFG2002(sys.argv[1], 1.0e-15, 300.0)

sfg.PrintData()
#sfg.PlotData()
#plt.show()
