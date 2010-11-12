import sys
from scipy import *
from scipy.optimize import leastsq, fmin_slsqp
import numpy

import matplotlib.pyplot as plt
import PlotUtility

from ColumnDataFile import ColumnDataFile as CDF

from StatisticMachine import StatisticMachine as SM

class LineShape:
	sm = SM()
	voigt_func = sm.voigt
	triple_voigt_func = sm.triple_voigt_fit

	lorentzian_func = sm.lorentzian
	triple_lorentzian_func = sm.triple_lorentzian_fit

	triple_funcs = {voigt_func:triple_voigt_func, lorentzian_func:triple_lorentzian_func}

	def __init__(self, file):
		self.cdf = CDF(file)
		self.x = self.cdf[0]
		self.data = self.cdf[1]

	def SetFunction(self,func):
		self.func = func
		self.triple_func = LineShape.triple_funcs[func]
	def SetParameters(self,p):
		self.p = p

	# smoothing is done with a moving average - gaussian moving filter
	def SmoothLine(self):
		self.smoothed = self.smoothListGaussian(self.data)
		min = self.x[0]
		max = self.x[-1]
		dx = (max-min)/len(self.smoothed)
		self.smoothed_x = [float(i)*dx+min for i in range(len(self.smoothed))]

	def FitLine(self):
		p0 = array (self.p)
		residual = LineShape.sm.FittingFunction(LineShape.sm.residuals,function=self.triple_func)
		plsq = leastsq(residual, p0, args=(self.data,self.x), maxfev=10000)
		print len(plsq)
		print plsq[0]

		self.plsq = plsq[0]
		self.fit = self.triple_func(self.x,self.plsq)
	
	def PlotLine(self,axs,lbl='Raw data'):
		axs.plot(self.x, self.data, linestyle=':', linewidth=2.0, label=lbl)
	def PlotSmoothLine(self,axs,lbl='Smooth data'):
		axs.plot(self.smoothed_x, self.smoothed, linestyle='-', color='r', linewidth=4.0, label=lbl)
	def PlotFit(self,axs,lbl=''):
		axs.plot(self.x, self.fit, linestyle='-', color='r', linewidth=3.5, label=lbl)

	def PlotComponentPeaks(self,axs):
		colors = ['orange', 'yellow', 'red', 'green']
		for i in range(3):
			# now calculate each of the 3 convolution peaks individually to show them off
			peak = self.func(self.x, self.plsq[i::3])
			axs.plot(self.x, peak, color='k', linewidth=1.0)
			axs.fill_between(self.x, 0, peak, color=colors[i], alpha=0.4)


	def smoothListGaussian(self, list, degree=20):  
		window=degree*2-1  
		weight=numpy.array([1.0]*window)  
		weightGauss=[]  
		for i in range(window):  
			i=i-degree+1  
	 		frac=i/float(window)  
	 		gauss=1/(numpy.exp((4*(frac))**2))  
	 		weightGauss.append(gauss)  
		weight=numpy.array(weightGauss)*weight  
		smoothed=[0.0]*(len(list)-window)  
		for i in range(len(smoothed)):  
			smoothed[i]=sum(numpy.array(list[i:i+window])*weight)/sum(weight)  
	
		return smoothed  

	def PrintParams(self):
		for i in range(len(self.plsq)/3):
			print ":\n\t%6.2f\n\t%6.2f\n\t%6.2f\n" % tuple(self.plsq[3*i:3*i+3])
		return





# get the data from the file

# peak parameters - first guesses for the fitting
triple_voigt_p = [3733, 3446, 3250,	# peak centers
								 5.0, 12.0, 12.0,
								 70.0, 100.0, 100.0,	# gaussian widths
								 2.0, 0.5, 2.0]			# amplitudes

triple_lorentzian_p = [3750, 3450, 3250,	# peak centers
											20.0, 20.0, 20.0,	# lorentzian widths
											1.0e-15, 2.0e-16, 1.0e-18]			# amplitudes

# plot the data
fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)

l = LineShape(sys.argv[1])
l.SetFunction(LineShape.lorentzian_func)
#l.SetFunction(LineShape.voigt_func)
#l.SetParameters(triple_voigt_p)
l.SetParameters(triple_lorentzian_p)
l.SmoothLine()
l.FitLine()
l.PlotLine(axs)
l.PlotSmoothLine(axs)
#l.PlotFit(axs)
#l.PlotComponentPeaks(axs)


PlotUtility.ShowLegend(axs)
plt.show()
