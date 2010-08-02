from ColumnDataFile import ColumnDataFile as CDF
from DataSetAnalyzer import DataSet
import sys
import matplotlib.pyplot as plt

def BondColors(bondlengths):
	colors = []
	for bond in bondlengths:
		val = 0
		for b in bond:
			if b < 2.46:
				val = val + 1
		colors.append(val)
	return colors
 
def PlotColors(ax,colors):
	x_last = 0
	x_current = 0
	color = 0
	for c in colors:
		x_current = x_current + 1
		if c != color:
			fc = 'g'
			if color == 2:
				fc = 'y'
			if color == 3:
				fc = 'r'
			if color != 0:
				ax.axvspan(x_last,x_current,alpha=0.3,facecolor=fc)
			color = c
			x_last = x_current
	

datacols = [3,5,7,11,13,15]
data = CDF(sys.argv[1],datacols)
dataset = [DataSet(data[col]) for col in datacols]

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)

for i in range(3):
  	ax1.plot(dataset[i].Data())
	ax2.plot(dataset[i+3].Data())

# do some cool coloring on the graph to show when 1 (green), 2 (yellow) and 3 (red) h-bonds are formed to a given oxygen
o1_bonds = zip(data[datacols[0]],data[datacols[1]],data[datacols[2]])
o2_bonds = zip(data[datacols[3]],data[datacols[4]],data[datacols[5]])

o1_colors = BondColors(o1_bonds)
o2_colors = BondColors(o2_bonds)

PlotColors(ax1,o1_colors)
PlotColors(ax2,o2_colors)

plt.show()
