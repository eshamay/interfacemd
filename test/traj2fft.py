#!/usr/bin/python

import sys
from TCFFile import TCFFile

import matplotlib.pyplot as plt

	
def PlotData(tcf):
	
  n = len(tcf)
  for i in range(n):
	chi_ax = fig.add_subplot(3,2,i+1)
#parts_ax = fig.add_subplot(n,1,2)

	chi_ax.plot(tcf[i].freq, tcf[i].chi, linewidth=2, label=tcf[i].filename)
#parts_ax.plot(tcf.freq, tcf.re, 'r-', linewidth=2)
#parts_ax.plot(tcf.freq, tcf.im, 'b-', linewidth=2)

#for ax in fig.get_axes():
	chi_ax.set_xlim(700.0, 3000.0)
	  #ax.set_ylim(0.0, 400.0)

	chi_ax.set_ylim(0.0,1.1)
#parts_ax.set_ylim(-800.0,800.0)
  	plt.legend()
	
def SetLegend(names):
	# set some legend properties.  All the code below is optional.  The
	# defaults are usually sensible but if you need more control, this
	# shows you how
	#
  	print names
	leg = plt.legend()
	# the matplotlib.patches.Rectangle instance surrounding the legend
	frame  = leg.get_frame()
	frame.set_facecolor('0.80')    # set the frame face color to light gray

	# matplotlib.text.Text instances
	for t in leg.get_texts():
		t.set_fontsize('medium')    # the legend text fontsize

	# matplotlib.lines.Line2D instances
	for l in leg.get_lines():
		l.set_linewidth(1.5)  # the legend line width



temps = [300.0,300.0,300.0,273.0,273.0]
tcf = map(lambda i,j: TCFFile(i,j), sys.argv[1:],temps)

# Set up the plot parameters (labels, size, limits, etc)
fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
PlotData(tcf)
#SetLegend([i.filename for i in tcf])
plt.show()
