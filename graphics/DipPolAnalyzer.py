'''
A class to take in a list of dipoles and a list of polarizabilities for a number of timesteps. The set of values is then manipulated to create various TCFs for certain sfg polarizations using the p,q,r (lab-frame coordinates) schema
'''
from numpy import *
from Utility import *

class DipPolAnalyzer:
	def __init__(self,rho,alpha=None):
	  	#print rho[0][0]
		self.rho = [array(r) for r in rho]
		#print self.rho[:10]
		if alpha != None:
			self.alpha = [apply(mat,[group(a,3)]) for a in alpha]

	def Rho(self):
		return self.rho
	def Alpha(self):
		return self.alpha

	def IR_TCF (self):
		ret = []
		try: 
			ret = [dot(i,self.rho[0]) for i in self.rho]
		except ValueError:
			print self.rho[0]

		return ret

	def SFG_TCF(self,s1,s2,p):
		alpha = [(a[s1,s1] + a[s2,s2])/2.0 for a in self.alpha]
		#alpha = [a[s1,s1] for a in self.alpha]
		rho = self.rho[0][p]
		return [a*rho for a in alpha]
		#return [dot(a,self.rho[0][r]) for a in self.alpha]

