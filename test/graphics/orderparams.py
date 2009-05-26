#!/Library/Frameworks/Python.framework/Versions/Current/bin/python

import sys
from OrderParameters import *

o = OrderParameters(sys.argv[1])
o.PlotData()
