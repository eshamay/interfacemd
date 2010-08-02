import matplotlib.pyplot as plt

class BondLengthFile:
  	def __init__(self,file):
		self.filename = file
		self.ParseFileContents(file)

  	def ParseFileContents(self,file):
		fp = open(file)
		self.data = [line.strip().split() for line in fp.readlines()]
		self.data = [map(lambda x: float(x),i) for i in self.data]

		self.x = [i[0] for i in self.data]
		self.y1 = [i[1] for i in self.data]
		self.y2 = [i[2] for i in self.data]

		fp.close()

	def PlotData(self,ax):

		ax.plot(x, y1, linewidth=1)
		ax.plot(x, y2, linewidth=1)

		#ax.set_xlim(700.0, 3000.0)
		#ax.set_ylim(0.0,1.1)
	def X(self):
	  	return self.x
	def Y1(self):
	  	return self.y1
	def Y2(self):
	  	return self.y2

