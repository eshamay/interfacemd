#!/usr/bin/python
import sys
from TCFFile import TCFFile

import matplotlib.pyplot as plt

class MoritaSFG2002:

	def __init__(self, file):
		## first open up the bond-length trajectory file
	  	self.tcf = TCFFile (file, 298.0)

	def PlotData(self):

		# Set up the plot parameters (labels, size, limits, etc)
		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

		# the extents of the x-range
		xmin = 2000.0
		xmax = 4000.0

		axs = fig.add_subplot(2,1,1)
		axs2 = fig.add_subplot(2,1,2)

  		for f in self.tcf.fourier:
			chi = [(c*c.conjugate()).real for c in f]
			chi_max = max(chi)
  			chi = [c/chi_max for c in chi]
			axs.plot(self.tcf.freq, chi, label='chi^2')
  			axs2.plot(self.tcf.freq, [c.real for c in f], label='real')
  			axs2.plot(self.tcf.freq, [c.imag for c in f], label='imaginary')

		axs.set_xlim(1000,4000)
		#axs.set_ylim(0.0, 0.01)
		axs2.set_xlim(1000,4000)




sfg = MoritaSFG2002(sys.argv[1])
sfg.PlotData()
plt.show()
