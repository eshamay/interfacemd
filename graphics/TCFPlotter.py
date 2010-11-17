from TCFSFGAnalyzer import TCFSFGAnalyzer 
from ColumnDataFile import ColumnDataFile as CDF
import sys

import matplotlib.pyplot as plt

cdf = CDF(sys.argv[1])
tcf = TCFSFGAnalyzer(cdf[0])
tcf.CalcSFGChi(1.0e-15, 298.0)

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
ax = fig.add_subplot(1,1,1)

ax.plot (tcf.Freq(), tcf.ChiSquared())
ax.set_xlim(1000,4000)
#ax.set_ylim(0.0,10e42)

plt.show()
