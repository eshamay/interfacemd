#!/usr/bin/python
import sys
from RDFFile import RDFFile
from RDF2DFile import RDF2DFile
from pylab import *

rdf = RDFFile (sys.argv[1])
rdf2d = RDF2DFile (sys.argv[2])

fig = figure(1)
ax = fig.add_subplot(1,1,1)

rdf.Plot(rdf.data, ax)
rdf2d.Plot1D(rdf2d.data, ax, [68.0, 70.0, 72.0])	# the hilite line should be red

axis([1.5, 11.0, 0.0, 6.0])
plt.show()

