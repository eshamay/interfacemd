'''
A class to take in a list of dipoles and a list of polarizabilities for a number of timesteps. The set of values is then manipulated to create various TCFs for certain sfg polarizations using the p,q,r (lab-frame coordinates) schema
'''
from numpy import *
from Utility import *

class DipPolAnalyzer:
	def __init__(self,rho,alpha):
	  	#print rho[0][0]
		self.rho = [array(r) for r in rho]
		#print self.rho[:10]
		self.alpha = [apply(mat,[group(a,3)]) for a in alpha]

	def IR_TCF (self):
		return [dot(i,self.rho[0]) for i in self.rho]

	def SFG_TCF(self,p,q,r):
		return [dot(a[p,q],self.rho[0][r]) for a in self.alpha]

