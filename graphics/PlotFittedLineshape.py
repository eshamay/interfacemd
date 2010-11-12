import sys
import matplotlib.pyplot as plt
import PlotUtility

from LineShape import LineShape

from scipy import *

# plot the data
fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.set_ylabel(r'# of SO$_2$ adsorbed', size='xx-large')
axs.set_xlabel('Time (ps)', size='xx-large')

l = LineShape(sys.argv[1])

l.PlotLine(axs)

l.SmoothLine()
l.PlotSmoothLine(axs)

l.SetParameters([30.0, 1.0e-3, 1.0e-4, -25.0])
l.SetFunction(LineShape.equil_rxn_func)
l.FitLine()
l.PlotFit(axs,clr='k',lbl='Equilibrium rxn',stl='-')

l.SetParameters([30.0, 1.0e-3, -25.0])
l.SetFunction(LineShape.forward_rxn_func)
l.FitLine()
l.PlotFit(axs,clr='r',lbl='Forward rxn',stl='--')

PlotUtility.ShowLegend(axs)
labels = axs.get_xticklabels() + axs.get_yticklabels()
for label in labels:
	label.set_size('xx-large')

plt.show()
