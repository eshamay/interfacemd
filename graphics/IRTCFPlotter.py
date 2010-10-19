#!/usr/bin/python
# assumes timestep of 1.0 fs and temp of 300K

import sys

from ColumnDataFile import ColumnDataFile
from DipPolAnalyzer import DipPolAnalyzer
from TCFSFGAnalyzer import TCFSFGAnalyzer

from PlotUtility import *
import matplotlib.pyplot as plt

class TCFIRPlotter:
	def __init__(self,file,dt=1.00e-15, temp=300.0):
		self.datafile = ColumnDataFile(file)
		self.dipoles = apply (zip, [self.datafile[i] for i in range(3)])
		self.dpa = DipPolAnalyzer(self.dipoles, None)

		self.dt = dt
		self.temp = temp
  	
	def PlotData(self):
		# Set up the plot parameters (labels, size, limits, etc)
		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
	  	
		# the extents of the x-range
		xmin = 500.0
		xmax = 4000.0

		# plotting out IR spectra
		ir = TCFSFGAnalyzer (self.dpa.IR_TCF())
		ir.CalcSFGChi(self.dt,self.temp)
		axs = fig.add_subplot(1,1,1)
		axs.plot(ir.Freq(), ir.ChiSquared(), label='IR')
		axs.set_xlim(xmin,xmax)
		ShowLegend(axs)



'''
class DipoleVectorIRPlotter:
	def __init__(self,file):
	  	self.datafile = ColumnDataFile(file)
		self.dipoles = apply (zip, [self.datafile[i] for i in range(3)])
		self.dpa = DipPolAnalyzer (self.dipoles, [])
		self.ir = TCFSFGAnalyzer(self.dpa.IR_TCF())
		self.ir.CalcSFGChi(0.75e-15,300.0)

	def PlotData(self):
		# Set up the plot parameters (labels, size, limits, etc)
		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
	  	
		# plotting out IR spectra
		axs = fig.add_subplot(1,1,1)
		axs.plot(self.ir.Freq()[1000:], self.ir.ChiSquared()[1000:], label='IR')
		axs.set_xlim(1000,4200)
'''

ir = TCFIRPlotter(sys.argv[1])
ir.PlotData()
plt.show()
