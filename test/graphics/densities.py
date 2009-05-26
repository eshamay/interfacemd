#!/usr/bin/python
import sys
from DensityProfiler import DensityProfiler

files = []
for file in sys.argv:
	files.append(file)
files = files[1:]

average = False
if not files[0].find('avg') == -1:
	average = True

d = DensityProfiler(files,avg=False)
d.PlotData()
#p0 = d.fits[0]['p_h2o']
#print "Lower interface is : % 8.3f\nUpper interface is : % 8.3f\n" % (p0[2], p0[4])
