#!/usr/bin/python
import sys
from ColumnDataFile import ColumnDataFile
from TCFSFGAnalyzer import TCFSFGAnalyzer

import matplotlib.pyplot as plt

class MoritaSFG2002:

	def __init__(self, file, dt=0.75e-15):
		## first open up the bond-length trajectory file
	  	self.tcf = ColumnDataFile(file,1)
  		self.sfg = [TCFSFGAnalyzer (self.tcf[key]) for key in self.tcf.keys()]
		map (lambda x: x.CalcFFT(dt), self.sfg)

	def PlotData(self,plotlist=[0]):

		# Set up the plot parameters (labels, size, limits, etc)
		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		# the extents of the x-range
		xmin = 2000.0
		xmax = 4000.0

		i = 1
		for p in plotlist:
			axs = fig.add_subplot(len(plotlist),1,i)
  			i = i + 1
			axs.plot(self.sfg[p].Freq()[1000:], [(x*x.conjugate()).real for x in self.sfg[p].FFTNormSquared()][1000:], label='chi^2')
			axs.set_xlim(1000,4000)


sfg = MoritaSFG2002(sys.argv[1])

sfg.PlotData([0,1,2,3,4,5,6,7,8,9,10,11])
plt.show()
