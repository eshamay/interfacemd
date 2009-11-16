#!/usr/bin/python
import sys
from DensityProfiler import DensityProfiler as dp

d = dp(sys.argv[1:])

d.PlotData(cols=[1])

