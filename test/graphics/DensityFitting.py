	def FitWater(self,data):
		x = data['position']
		h2o = data['h2o']

		# initial guess values
		p1 = 0.0
		p2 = 1.0
		p3 = 1.0
		p4 = 0.0
		x0 = 25.0
		d0 = 5.0
		x1 = 70.0
		d1 = 5.0

		if self.avg:
			p0 = array ([p1,p2,x0,d0])
			fit_func = self.sm.tanh_fit
		else:
			p0 = array ([p1,p2,p3,p4,x0,d0,x1,d1])
			fit_func = self.sm.double_tanh_fit

		residual=self.sm.FittingFunction(self.sm.residuals,function=fit_func)
		plsq = leastsq(residual, p0, args=(h2o,x), maxfev=100000)

		fit = fit_func(x,plsq[0])

		return (fit, p0)


