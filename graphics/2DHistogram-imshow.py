#!/usr/bin/python
import sys
from scipy import *
import scipy.io.array_import

# plotting libs
import matplotlib.pyplot as plt 
from pylab import *


# Load up multiple files and parse the matrix data from them
data = []
for arg in sys.argv[1:]:
	file = scipy.io.array_import.read_array(arg)
	datum = []
	for i in range(len(file[0])):
		datum.append(file[:,i])
	data.append(datum)

# create a figure
fig = figure(1)

cmBase = cm.jet   # choose a colormap to use for the plot

# This processes the data and adds the colored 2D histogram to the figure
'''
for datum in data:
	im = plt.imshow(datum, aspect='auto', extent=(-20.0, 150.0, -1.0,1.0), interpolation='bilinear')    
'''
# Create the layout of the plots
titles = ['Bisector','Normal',r'OH$_A$',r'OH$_B$']
for i in range(len(data)):
	subplot(2,2,i+1)
	title(titles[i], fontsize=28)

	interface = 0.0
	# top row
	if i <= 1:
		interface = 65.744
	else:
		interface = 67.563


	im = plt.imshow(data[i], cmap=cm.jet, aspect='auto', extent=(-20.0-interface, 150.0-interface, 1.0, -1.0), interpolation='bilinear')    

  	# Left column
	# Turn on the X-Y grid
	if i%2 == 0:
		ylabel(r'$\cos\theta$', fontsize=28)
		yticks(fontsize=20)
	else:
		ylabel(r'$\cos\phi$', fontsize=28)
		yticks(yticks()[0], ())

  	# bottom row
	if i > 1:
		xlabel('Slab Position', fontsize=28)
		xticks(fontsize=20)
	else:
		xticks(yticks()[0], ())


	grid(True)

for a in fig.axes:
	a.set_xlim(-3.0, 7.0)

# A little magic to create a colorbar - legend for the plot(s)
cax = axes([0.95, 0.33, 0.025, 0.33])
cb = plt.colorbar(im, cax=cax) # grab the Colorbar instance
for t in cb.ax.get_yticklabels():
	t.set_fontsize(20)
cb.ax.set_yticklabels(())	

# Pass go - collect $200
draw()
show()
