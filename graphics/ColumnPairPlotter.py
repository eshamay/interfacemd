import sys
import matplotlib.pyplot as plt
from ColumnPairFile import ColumnPairFile as CPF
import PlotUtility

from matplotlib import collections, axes, transforms
from matplotlib.colors import colorConverter
import numpy as N

def PlotOffsetColumnPairs (files, column=None):

	cpf = [CPF(f) for f in files]
	fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
	#axs = [fig.add_subplot(numPlots,1,i+1) for i in range(numPlots)]
	axs = fig.add_subplot(1,1,1)
	
	segs = []
	for f in cpf:
		#f.PlotData(axs, column)
		p = f.ColumnPair(column)
		print f.filename
		line = zip(p[0],p[1])
		segs.append(line)

	col = collections.LineCollection(segs, offsets=(0.0,1.0))
	axs.add_collection(col, autolim=True)
	axs.autoscale_view()

	plt.show()
