#!/usr/bin/python
import sys
from SimpleDensityProfiler import DensityProfiler as dp

d = dp(sys.argv[1])

d.PlotData()
