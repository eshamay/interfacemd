import sys
import matplotlib.pyplot as plt
from ColumnPairFile import ColumnPairFile as CPF


numPlots = int(sys.argv[1])
cpf = [CPF(f) for f in sys.argv[2:]]
fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = [fig.add_subplot(numPlots,1,i+1) for i in range(numPlots)]

for f in cpf:
	for a in range(numPlots):
		f.PlotData(axs[a], a)

for a in axs:
	a.set_ylim(0,1.2)

plt.show()
