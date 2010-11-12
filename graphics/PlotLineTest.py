from StatisticMachine import StatisticMachine as SM
import matplotlib.pyplot as plt
import PlotUtility

from scipy import *
import PlotUtility

sm = SM()
# plot the data
fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)

x = array(range(0,1500000,250))
fcold = sm.rate_rxn(x,[  6.46175689e+01,   3.03049349e-06,   2.69925293e-06,  -2.31903106e+01])
fwarm = sm.rate_rxn(x,[  7.43871131e+01,   2.90845091e-06,   2.24327538e-06,  -2.50737245e+01])

axs.plot(x,fwarm, label='warm')
axs.plot(x,fcold, label='cold')
PlotUtility.ShowLegend(axs)
plt.show()
