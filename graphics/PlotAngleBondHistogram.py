import sys
import Utility
import PlotUtility
from ColumnPairFile import ColumnPairFile as CPF

import matplotlib.pyplot as plt

file = sys.argv[1]
cpf = CPF(file)
fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
cpf.PlotData(ax1,0)
cpf.PlotData(ax2,1)
ax1.set_xlabel('Bondlength', size='xx-large')
ax2.set_xlabel('Angle', size='xx-large')
ax2.set_xlim(0,180)
plt.show()
