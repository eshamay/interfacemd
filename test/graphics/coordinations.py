#!/usr/bin/python
import sys
from Coordinator import Coordinator as coord

files = []
for file in sys.argv[1:]:
	files.append(file)

c = coord(files,norm=True)
c.PlotData()
