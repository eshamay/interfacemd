import Utility
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
from matplotlib import collections

def ScaleToMax(data):
	data_max = max(data)
	data = map(lambda x: x/data_max, data)
	return data

files = [i for i in Utility.SearchDirectoryTree('.','cp2k-morita2002.xvg')]
for f in files:
	print f

cdfiles = [CDF(i) for i in files]

titles = ['IR', 'SFG X', 'SFG Y', 'SFG Z']
for i in range(1,5):
	# get the collection of lines assembled
	line = [zip(cd[0], ScaleToMax(cd[i])) for cd in cdfiles]

	col = collections.LineCollection(line, offsets=(0.0,1.0))
	fig = plt.figure(num=i, facecolor='w', edgecolor='w', frameon=True)
	axs = fig.add_subplot(1,1,1, title=titles[i-1])
	axs.add_collection(col, autolim=True)
	axs.autoscale_view()
	axs.set_xlim(0,5000)

plt.show()
