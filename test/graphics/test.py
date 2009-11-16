import sys
from scipy import *
import scipy.io.array_import

# plotting libs
import matplotlib.pyplot as plt 
from pylab import *


# Load up three different files
file1 = scipy.io.array_import.read_array(sys.argv[1])
file2 = scipy.io.array_import.read_array(sys.argv[2])
data1 = []
data2 = []

for i in range(len(file1[0])):
		data1.append(file1[:,i])
for i in range(len(file2[0])):
		data2.append(file2[:,i])
x = data1[0]
y =	data1[1]

# create a figure
fig = plt.figure(1, facecolor='white', figsize=(11,11))
# Create the triptych layout of the three histograms
ax = fig.add_subplot(1,1,1)
plt.title('Molecular Axis Orientation Distribution', fontsize=28)
ax.plot (x, y, 'r-', linewidth=4, label='Decane')

x = data2[0]
y =	data2[1]
ax.plot (x, y, 'b-', linewidth=4, label='Perfluorodecane')

plt.legend()

# zoom in on the region of interest
#plt.ylim(2.5, 32.5)
#plt.xlim(-1.5, 30.0)
# Turn on the X-Y grid
yticks(fontsize=18)
xlabel(r'$\cos(\theta)$', fontsize=28)		
ylabel(r'$\rho$', fontsize=28)
xticks(fontsize=18)

grid(True)

plt.show()
