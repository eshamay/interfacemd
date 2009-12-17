import sys
from scipy import *
import scipy.io.array_import

# plotting libs
import matplotlib.pyplot as plt 
from pylab import *


# Load up three different files
data = []
for arg in sys.argv[1:]:
	file = scipy.io.array_import.read_array(arg)
	datum = []
	for i in range(len(file[0])):
		datum.append(file[:,i])
	data.append(datum)

# create a figure
fig = figure(1, facecolor='white')

cmBase = cm.jet   # choose a colormap to use for the plot

# This processes the data and adds the colored 2D histogram to the figure
'''
for datum in data:
im = plt.imshow(datum, aspect='auto', extent=(-20.0, 150.0, -1.0,1.0), interpolation='bilinear')    
'''
# Create the triptych layout of the three histograms
titles = ['Bisector','Normal',r'OH$_A$',r'OH$_B$']
for i in range(len(data)):
	subplot(2,2,i+1)
	title(titles[i], fontsize=28)
	im = plt.imshow(data[i], cmap=cm.jet, aspect='auto', extent=(-1.0, 1.0, 150.0, -20.0), interpolation='bilinear')    

	# zoom in on the region of interest
	plt.ylim(68.0, 58.0)
	# Turn on the X-Y grid
	if i%2 == 0:
		yticks(fontsize=20)
		ylabel('Slab Position', fontsize=28)
	else:
		yticks(yticks()[0], ())
	if i > 1:
		xlabel(r'$\cos(\theta)$', fontsize=28)		
		xticks(fontsize=20)
		
	grid(True)

# A little magic to create a colorbar - legend for the plot(s)
cax = axes([0.95, 0.33, 0.025, 0.33])
cb = plt.colorbar(im, cax=cax) # grab the Colorbar instance
for t in cb.ax.get_yticklabels():
	t.set_fontsize(20)
cb.ax.set_yticklabels(())	

# Pass go - collect $200
draw()
show()
