from scipy import *

class StatisticMachine:

	def __init__(self):
		pass

	# calculate the residual error between the two linesets (experimental and simulated)
	# the exp & sim params are tuples of (x,y) data
	def residual(self,exp,sim):

		res = []
		for i in range(len(exp[0])):
			res.append(exp[1][i] - sim[1][i])

		return res

	# sum of squares of the regression
	def SSR(self,exp,sim):

		SSR = 0.0
		mean = self.mean(exp[1])

		for d in sim[1]:
			SSR = SSR + (d - mean)**2
		return SSR

	# sum of squares due to error
	def SSE(self,exp,sim):

		res = self.residual(exp,sim)

		SSE = 0.0
		for r in res:
			SSE = SSE + r**2

		return SSE

	# sum of squares about the mean
	def SST(self,exp,sim):

		mean = self.mean(exp[1])
		SST = 0.0

		for d in exp[1]:
			SST = SST + (d - mean)**2

		return SST

	# R-squared of a fit
	def r_square(self,exp,sim):

		SSE = self.SSE(exp,sim)
		SST = self.SST(exp,sim)

		return 1 - SSE/SST

	# calculate the mean of a dataset
	def mean(self,data):
		y_bar = 0.0
		for y in data:
			y_bar = y_bar + y
		y_bar = y_bar / len(data)

		return y_bar


	def mean_square(self,data):
		y_bar = 0.0
		for y in data:
			y_bar = y_bar + y**2
		y_bar = y_bar / len(data)

		return y_bar

	# fitting functions for trigonometric and gaussian lines

	# parameters are:
	# 	p = parameter list
	# 	y = y, x = x,
	# 	function = fitting function to be used
	def residuals(self,p, y, x, function):
		err = y - function(x,p)
		return err

	def FittingFunction(self, f, **kwargs):
		return lambda *largs: f(*largs,**kwargs)

	# 0,1 = rho1 & rho2
	# 2,3 = Z0, d
	def tanh_fit(self,x,p):
		return 0.5*(p[0]+p[1]) - 0.5*(p[0]-p[1])*tanh((x-p[2])/p[3])

	# 0 = density of side1
	# 1 = density of side2
	# 2,3 = density/rho sides 3 & 4
	# 4 = z0 location of first interface
	# 5 = d of 1st int.
	# 6,7 = z0 and d of 2nd interface
	# 8,9,10 = a,b,c of first ion
	# 11,12,13 = a,b,c of 2nd ion
		# a = height of gaussian
		# b = offset
		# c = width
	# used for unaveraged (double) interfaces
	def double_tanh_gaussian_fit(self,x,p):
		return 0.5*(p[0]+p[1]) - 0.5*(p[0]-p[1])*tanh((x-p[4])/p[5]) - 0.5*(p[2]-p[3])*tanh((x-p[6])/p[7]) + p[8]*exp(-((x-p[9])**2)/(2.0*p[10]**2)) + p[11]*exp(-((x-p[12])**2)/(2.0*p[13]**2))


	# 0 = density of side1
	# 1 = density of side2
	# 2,3 = Z0,d of the tanh
	# 4,5,6 = a,b,c of the gaussian
	def tanh_gaussian_fit(self,x,p):
		### here's your standard normal (gaussian) distribution
		return 0.5*(p[0]+p[1]) - 0.5*(p[0]-p[1])*tanh((x-p[2])/p[3]) + p[4]*exp(-(x-p[5])*(x-p[5])/2.0/p[6]/p[6])
		### this is a try at the lorentzian... which doesn't work!
		#return 0.5*(p[0]+p[1]) - 0.5*(p[0]-p[1])*tanh((x-p[2])/p[3]) + 1/pi*(p[4]/((x-p[5])**2 + p[4]**2))

	# 0,1 = rho1 & rho2 of side 1
	# 2,3 = rho1 & rho2 of side 2
	# 4,5 = Z0,d side 1
	# 6,7 = Z0,d side 2
	def double_tanh_fit(self,x,p):
		return 0.5*(p[0]+p[1]) - 0.5*(p[0]-p[1])*tanh((x-p[4])/p[5]) - 0.5*(p[2]-p[3])*tanh((x-p[6])/p[7])
