"""
Utilities
"""
from math import *
from PhysicalConstants import *

def wait():
    raw_input("Press a key...")


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step
def lin(x,x0,y0,x1,y1):
	return y0 + (x-x0)*(y1-y0)/(x1-x0)


class Point3D:
	def __init__(self,x=0,y=0,z=0):
		"""
		A 3D point  
		
		"""
		
		self.x = x  # creation 
		self.y = y
		self.z = z
		
	def X(self):
		"""
		x coordinate 
		"""
		return self.x

	def Y(self):
		"""
		Y coordinate 
		"""
		return self.y

	def Z(self):
		"""
		z coordinate 
		"""
		return self.z

	def XYZ(self):
		"""
		xyz coordinates
		"""
		return self.x,self.y,self.z

	def __str__(self):

		s= """
	      x = %7.2f  y = %7.2f  z = %7.2f  
			"""%(self.X(), self.Y(), self.Z())

		return s

class Point4D(Point3D):
	def __init__(self,x=0,y=0,z=0, t=0):
		"""
		A 4D point  
		
		"""
		Point3D.__init__(self,x,y,z)
		self.t = t

	def T(self):
		"""
		t coordinate 
		"""
		return self.t

	def XYZT(self):
		"""
		t coordinate 
		"""
		return self.x,self.y,self.z,self.t

	def __str__(self):

		s= """
	      x = %7.2f  y = %7.2f z = %7.2f  t = %7.2f  
			"""%(self.X(), self.Y(), self.Z(), self.T())

		return s 


if __name__ == '__main__':
	p = Point3D(10,10,-400)
	print p
	p = Point4D(10,10,-400,2)
	print p

