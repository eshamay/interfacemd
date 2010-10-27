import sys
import Utility
from ColumnPairFile import ColumnPairFile as CPF

from PlotUtility import *
import matplotlib.pyplot as plt

Na = 6.02e23
dr = 0.1
size = [40.0, 60.0, 40.0]

densities = CPF(sys.argv[1],True)

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

for d in range(len(densities)):
	size_scale = size[1]*size[2]*size[0]
	if d > 2:
		size_scale = size_scale / size[d%3]
	else:
		size_scale = size_scale / size[d]
	ax = fig.add_subplot(len(densities)/3,3,d+1)
	#scale = 0.33241	# 18.01 / Na / 1e-24 / xsize / ysize
	scale = float(sys.argv[2]) / dr / size_scale / .602
	densities.PlotData(ax,d,scale2=scale)
	ShowLegend(ax)

plt.show()


