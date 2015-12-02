import random as rnd
from LXe import *
from Geometry import *
from Util import *
from Centella.histoManager import *
from Centella.messenger import *


"""
Generates photons in a point located at coordinates (0,0,0) in a system
of reference in which box1 and box2 are located as in the description below.
The box coordinates are defined by 8 vertices.
Each vertex is an (x,y,z) point. 
Define the box with the following convention:
        
        v1 = (xl,yb,zb), where xl = lefmost x coordinate,
                               yb = bottom y coordinate
                               zb = back z coordinate
        x-------x (xr,yt)
        |       |
        |       |
        x-------x            zb ----- zf
        (xl,yb)

        v2 = (xl,yt,zb), where xl = lefmost x coordinate,
                               yt = top y coordinate
                          

        v3 = (xr,yb,zb), where xr = rightmost x coordinate
        v4 = (xr,yt,zb)
        v5 = (xl,yb,zf), where zf = front z coordinate
        v6 = (xl,yt,zf)
        v7 = (xl,yt,zf)
        v8 = (xr,yb,zf)                           


"""

BOX1 =[
[-12.8,-12.8,-100.],[-12.8,12.8,-100.],[12.8,-12.8,-100.],[12.8,12.8,-100.],
[-12.8,-12.8,-130.],[-12.8,12.8,-130.],[12.8,-12.8,-130.],[12.8,12.8,-130.]
]

BOX2 =[
[-12.8,-12.8,100.],[-12.8,12.8,100.],[12.8,-12.8,100.],[12.8,12.8,100.],
[-12.8,-12.8,130.],[-12.8,12.8,130.],[12.8,-12.8,130.],[12.8,12.8,130.]
]

class PhotonGenerator:
	def __init__(self,boxCoord1,boxCoord2,nevents,level):
		"""
		boxCoord1,boxCoord2 := coordiantes of the boxes
		nevents = number of events to generate
		level = debug level
		 
		
		"""
		self.box1 = Box(boxCoord1)
		self.box2 = Box(boxCoord2) 
		self.m = Messenger(level)
		self.hman =HistoManager() 
		self.nevents = nevents
		self.m.log(2, "Box1 --", self.box1)
		self.m.log(2, "Box2 ---", self.box2)

		self.lxe = LXe() #lxe properties
		self.m.log(2, "LXe ---", self.lxe)

		
		
	def GenerateMomentum(self):
		"""
		generate (px,py,pz) for photon 1
		"""
		

	
	def __str__(self):

		s= """
		Photon Generator:
		  box 1  = %s
		  box 2 = %s
		  number of events to generate = %d
	      
        
			"""%(self.box1,self.box2,self.nevents
				)

		return s



if __name__ == '__main__':
	pg = PhotonGenerator(BOX1,BOX2,100,4)
	print pg
	
	

