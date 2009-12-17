import sys
from scipy import *
import scipy.io.array_import

# plotting libs
import matplotlib.pyplot as plt 
from pylab import *

# load the file with the 2D histogram (matrix) data
file = scipy.io.array_import.read_array(sys.argv[1])
data = []
for i in range(len(file[0])):
		data.append(file[:,i])
x = data[0]
y =	data[1]

# create a figure
fig = plt.figure(1, facecolor='white', figsize=(11,11))

# Create the triptych layout of the three histograms
ax = fig.add_subplot(1,1,1)
plt.title('Interface Carbons - Decane', fontsize=28)
plt.suptitle('Location and Height of Topmost Carbon Atoms', fontsize=20)
ax.plot (x, y, 'k-', linewidth=4)

# zoom in on the region of interest
#plt.ylim(2.5, 32.5)
#plt.xlim(-1.5, 30.0)
# Turn on the X-Y grid
yticks(fontsize=18)
xlabel(r'Interface Position / $\AA$', fontsize=28)		
ylabel(r'Interface Position / $\AA$', fontsize=28)
xticks(fontsize=18)

grid(True)

plt.show()
