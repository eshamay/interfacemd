#!/usr/bin/python
import sys
from MoritaSFG import MoritaSFG

file = sys.argv[1]
d = MoritaSFG(file)
d.PlotData()
