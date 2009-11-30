#!/usr/bin/python

import sys
from ColumnDataFile import ColumnDataFile as CDF

import matplotlib.pyplot as plt

data = CDF(sys.argv[1])
x = data['Distance']

# Set up the plot parameters (labels, size, limits, etc)
self.fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

# some of the prelim stuff for our figure
ax = self.fig.add_subplot(1,1,1)
#ax.set_title(TITLE, size='x-large')

for key in data.keys():
	datum = data[key]
	ax.plot(x, datum, linewidth=1, label=key)
