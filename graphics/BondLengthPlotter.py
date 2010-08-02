#!/usr/bin/python
import sys
from BondLengthFile import BondLengthFile
from DataSetAnalyzer import DataSet

import matplotlib.pyplot as plt

blf = BondLengthFile(sys.argv[1])
bond_1 = DataSet(blf.y1)
bond_2 = DataSet(blf.y2)

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
length_ax = fig.add_subplot(2,1,1)
fft_ax = fig.add_subplot(2,1,2)

length_ax.plot(blf.X(), bond_1.data)
length_ax.plot(blf.X(), bond_2.data)

bond_1.CalcFFT(0.75e-15)
bond_2.CalcFFT(0.75e-15)
fft_ax.plot(bond_1.Freq()[20:], bond_1.FFTNormSquared()[20:])
fft_ax.plot(bond_2.Freq()[20:], bond_2.FFTNormSquared()[20:])

plt.xlim(500,4000)
plt.show()

# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
