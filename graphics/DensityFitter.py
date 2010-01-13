from scipy import *
from scipy.optimize import leastsq

from StatisticMachine import StatisticMachine as SM

class DensityFitter(SM):

	def __init__(self):
		pass

	def FitWater(self,x,data,gibbs_1=30.0,gibbs_2=70.0,avg=False):
		# initial guess values
		dens_low_1 = 0.0
		dens_high_1 = max(data)
		dens_high_2 = max(data)
		dens_low_2 = 0.0
		width_1 = 5.0
		width_2 = 5.0

		if avg:
			p0 = array ([dens_low_1,dens_high_1,gibbs_1,width_1])
			fit_func = self.tanh_fit
		else:
			p0 = array ([dens_low_1,dens_high_1,dens_high_2,dens_low_2,gibbs_1,width_1,gibbs_2,width_2])
			fit_func = self.double_tanh_fit

		residual = self.FittingFunction(self.residuals,function=fit_func)
		plsq = leastsq(residual, p0, args=(data,x), maxfev=100000)

		fit = fit_func(x,plsq[0])
		
		param_hash = {"dens_low_1":plsq[0][0], "dens_high_1":plsq[0][1], "dens_low_2":plsq[0][3], "dens_high_2":plsq[0][2], "gibbs_1":plsq[0][4], "width_1":plsq[0][5], "gibbs_2":plsq[0][6], "width_2":plsq[0][7]}

		print "Fitting the water density - found the following parameters:"
		print "Lower interface:"
		print "density (%6.3f) [%6.3f - %6.3f]    d = % 6.4f\n" % (plsq[0][4], plsq[0][0], plsq[0][1], plsq[0][5])

		print "Upper interface:"
		print "density (%6.3f) [%6.3f - %6.3f]    d = % 6.4f\n" % (plsq[0][6], plsq[0][3], plsq[0][2], plsq[0][7])

		print "Average of the densities (low-high) = %6.3f" % (sum([plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3]])/4.0)
		return (fit, param_hash)

	def FitLowerWater(self,x,data,gibbs=30.0):
		# initial guess values
		dens_low = 0.0
		dens_high = max(data)
		width = 5.0

		p0 = array ([dens_low,dens_high,gibbs,width])
		fit_func = self.tanh_fit

		residual = self.FittingFunction(self.residuals,function=fit_func)
		plsq = leastsq(residual, p0, args=(data,x), maxfev=100000)

		fit = fit_func(x,plsq[0])
		
		params = {"dens_low":plsq[0][0], "dens_high":plsq[0][1], "gibbs":plsq[0][2], "width":plsq[0][3]}

		print "Fitting the water density - found the following parameters:"
		print "Lower interface:"
		print "density (%6.3f) [%6.3f - %6.3f]    d = % 6.4f\n" % (plsq[0][2], plsq[0][0], plsq[0][1], plsq[0][3])

		return (fit, params)


