#!/usr/bin/python

import sys

from ColumnDataFile import ColumnDataFile
from TCFSFGAnalyzer import TCFSFGAnalyzer
from DipPolAnalyzer import DipPolAnalyzer
from PlotUtility import *
from numpy import dot
import matplotlib.pyplot as plt

class TCFIRPlotter:
	def __init__(self,file):
	  	self.datafile = ColumnDataFile(file)
		self.ir = TCFSFGAnalyzer(self.datafile[0])
		self.ir.CalcSFGChi(0.75e-15,300.0)
  	
	def PlotData(self):
		# Set up the plot parameters (labels, size, limits, etc)
		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
	  	
		# plotting out IR spectra
		axs = fig.add_subplot(1,1,1)
		axs.plot(self.ir.Freq()[1000:], self.ir.ChiSquared()[1000:], label='IR')
		axs.set_xlim(1000,4200)


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



if len(sys.argv) < 2:
  	print "Need both arguments <tcf/vec> <filename>"
	
else:
  	if sys.argv[1] == "vec":
		ir = DipoleVectorIRPlotter(sys.argv[2])
		ir.PlotData()
		plt.show()
	elif sys.argv[1] == "tcf":
		ir = TCFIRPlotter(sys.argv[2])
		ir.PlotData()
		plt.show()
	else:
	  	print "first argument is incorrect. Should be <tcf/vec>."

