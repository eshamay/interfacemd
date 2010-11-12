import sys
from scipy import *
from scipy.optimize import leastsq, fmin_slsqp
import numpy

from ColumnDataFile import ColumnDataFile as CDF

from StatisticMachine import StatisticMachine as SM

class LineShape:
	sm = SM()

	exp_inc_func = sm.exponential_decay_increasing
	lorentzian_func = sm.lorentzian
	voigt_func = sm.voigt
	logistic_func = sm.logistic
	steady_state_func = sm.steady_state
	equil_rxn_func = sm.equilibrium_rxn
	forward_rxn_func = sm.forward_rxn

	def __init__(self, file):
		self.cdf = CDF(file)
		self.x = self.cdf[0]
		self.data = self.cdf[1]

	def SetFunction(self,func):
		self.func = func

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
		residual = LineShape.sm.FittingFunction(LineShape.sm.residuals,function=self.func)

		plsq = leastsq(residual, p0, args=(self.data,self.x), maxfev=10000)

		self.plsq = plsq[0]
		print self.plsq
		self.fit = self.func(array(self.x),self.plsq)
	
	def PlotLine(self,axs,lbl='Raw data'):
		axs.plot(self.x, self.data, linestyle=':', linewidth=0.4, label=lbl)

	def PlotSmoothLine(self,axs,lbl='Smooth data'):
		axs.plot(self.smoothed_x, self.smoothed, linestyle='--', color='g', linewidth=2.0, label=lbl)

	def PlotFit(self,axs,lbl='',clr='r',stl='-'):
		axs.plot(self.x, self.fit, linestyle=stl, color=clr, linewidth=3.0, label=lbl)

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

