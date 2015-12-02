import random as rnd
from LXe import *
from Geometry import *
from Util import *
from Centella.histoManager import *
from Centella.treeManager import *
from Centella.messenger import *
from system_of_units import *
from NumpyAlgebra import *
import array 


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
[-12.8,-12.8,-150.],[-12.8,12.8,-150.],[12.8,-12.8,-150.],[12.8,12.8,-150.]
]

BOX2 =[
[-12.8,-12.8,100.],[-12.8,12.8,100.],[12.8,-12.8,100.],[12.8,12.8,100.],
[-12.8,-12.8,150.],[-12.8,12.8,150.],[12.8,-12.8,150.],[12.8,12.8,150.]
]
EPHOT = 511*keV
STEP = 1*mm

class PhotonGenerator:
	def __init__(self,box1,box2,level):
		"""
		box1, box2: instances of boxes
		level = debug level
		 
		"""
		self.m = Messenger(level)
		self.box1 = box1
		self.box2 = box2 
		self.xmin = self.box1.xmin #assume that box1 defines fiducial
		self.xmax = self.box1.xmax #assume that box1 defines fiducial
		self.dx = self.box1.x/2.
		self.ymin = self.box1.ymin #assume that box1 defines fiducial
		self.ymax = self.box1.ymax #assume that box1 defines fiducial
		self.dy = self.box1.y/2.
		self.dz = abs(self.box1.zmin) 
		
		self.m.log(2, "Box1 --", self.box1)
		self.m.log(2, "Box2 ---", self.box2)

		self.lxe = LXe() #lxe properties
		self.m.log(2, "LXe ---", self.lxe)

		xPE = self.lxe.PhotoelectricCrossSection(EPHOT)
		xTot = self.lxe.TotalCrossSection(EPHOT)
		ProbPE = xPE/xTot

		self.ProbStep = self.lxe.Efficiency(EPHOT,STEP)*ProbPE

		self.m.log(2, "Photoelectric probability for 511 keV = %7.2f"%(
			ProbPE))
		self.m.log(2, "Probability PE per mm at 511 keV = %7.2g"%(
			self.ProbStep))

	def GenerateEvent(self):
		"""
		Generates:
		gamma1 and gamma2 with normalized momentum P1 and P2
		by convention, P1 has negative Pz and P2 positive Pz
		"""
		P2 = self.generateMomentum(lvl=3)
		P1 = -1.*P2

		self.m.log(3," P1 =%s, P2 =%s"%(P1,P2))
			
		self.m.log(3," |P2| =%7.2f, |P1| =%7.2f, P1*P2 =%7.2f "%(
			P2.Norm(),P1.Norm(),P1*P2))

		#extrapolate P1 to box1 (z <0)
		z1 = self.box1.zmin
		x1,y1 = extrapToZ(z1,(0,0,0),(P1[0],P1[1],P1[2]))

		#extrapolate P2 to box2 (z >0)
		z2 = self.box2.zmin
		x2,y2 = extrapToZ(z2,(0,0,0),(P2[0],P2[1],P2[2]))

		
		self.m.log(3, 'x1 =%7.2f mm, y1 =%7.2f mm, z1 =%7.2f mm '%(
					x1/mm,y1/mm,z1/mm))

		self.m.log(3, 'x2 =%7.2f mm, y2 =%7.2f mm, z2 =%7.2f mm '%(
					x2/mm,y2/mm,z2/mm))

		#safety check

		if self.box1.Active((x1,y1,z1)) == False:
			self.m.log(0, '**generated point in box 1 out of box 1*** ')
			self.m.log(0, ' x1 =%7.2f mm, y1 =%7.2f mm, z1 =%7.2f mm '%(
					x1/mm,y1/mm,z1/mm))
			return False

		if self.box2.Active((x2,y2,z2)) == False:
			self.m.log(0, '**generated point in box 2 out of box 2*** ')
			self.m.log(0, ' x2 =%7.2f mm, y2 =%7.2f mm, z2 =%7.2f mm '%(
					x2/mm,y2/mm,z2/mm))
			return False

		#compute path in box1 and in box2

		path1 = pathInBox((x1,y1,z1),(P1[0],P1[1],P1[2]),
			self.box1)

		path2 = pathInBox((x2,y2,z2),(P2[0],P2[1],P2[2]),
			self.box2)
				

		self.m.log(3,'Prob of int for g1 (511 keV), path =%7.2f mm) = %7.2f'%(
				path1/mm,self.lxe.Efficiency(EPHOT,path1)))

		self.m.log(3,'Prob of int for g1 (511 keV), path =%7.2f mm) = %7.2f'%(
				path2/mm,self.lxe.Efficiency(EPHOT,path2)))

		#propagate photon 1:
		inter =0
		istep =0
		for step in drange(0, path1, STEP):
			self.m.log(4,'step =%7.2f mm '%(step/mm))
			#interacts?
			if rnd.uniform(0.,1.) <= self.ProbStep:
				inter =1
				istep = step
				break


		if inter == 0:
			return False  #fail in box1 no need to try in box2


		self.m.log(3,' Photon 1 interacts in step =%7.2f  '%(step))
		step1 = step

		#propagate photon 2:
		inter =0
		istep =0
		
		for step in drange(0, path2, STEP):
			self.m.log(4,'step =%7.2f mm '%(step/mm))
			#interacts?
			if rnd.uniform(0.,1.) <= self.ProbStep:
				inter =1
				istep = step
				break

		if inter == 0:
			return False  #fail in box2 no need to carry on

		self.m.log(3,' Photon 2 interacts in step =%7.2f  '%(step))
		step2 = step

		#finally, find the interaction point in each box.
		xi1,yi1,zi1 = propagateInBox((x1,y1,z1), 
			(P1[0],P1[1],P1[2]), step1)

		xi2,yi2,zi2 = propagateInBox((x2,y2,z2), 
			(P2[0],P2[1],P2[2]), step2) 

		self.m.log(2, 'xi1 =%7.2f mm, yi1 =%7.2f mm, zi1 =%7.2f mm '%(
					xi1/mm,yi1/mm,zi1/mm))

		self.m.log(2, 'xi2 =%7.2f mm, yi2 =%7.2f mm, zi2 =%7.2f mm '%(
					xi2/mm,yi2/mm,zi2/mm))

		self.X1 = [xi1,yi1,zi1]
		self.X2 = [xi2,yi2,zi2]
		return True
		

		
	def generateMomentum(self,lvl=0):
		"""
		generate (px,py,pz) for photon 1
		px and py are generated inside box solid angle
		tan(tx) = (x/2)/(z/2) = x/z
		where x is the length of the box and z the distance
		between the two boxes.
		px = (+-) E*tan(tx)
		where E = 511 keV and (+-) denotes a random sign
		"""

		
		px = self.generatePt(EPHOT,self.dx,self.dz)
		py = self.generatePt(EPHOT,self.dy,self.dz)
		pz = sqrt(EPHOT**2 - px**2 - py**2)

		self.m.log(lvl, " px = %7.2f py = %7.2f pz = %7.2f "%(
			px/EPHOT,py/EPHOT,pz/EPHOT))

		P = NPVector([px/EPHOT,py/EPHOT,pz/EPHOT]) 
		return P

	def generatePt(self,p,x,z):
		xf = x - 1*mm
		ctx = rnd.uniform(0.,xf/z)
		px = p*atan(ctx)
		if rnd.uniform(0.,1.) > 0.5: 
			px = -px
		
		return px

	def __str__(self):

		s= """
		Photon Generator:
		  box 1  = %s
		  box 2 = %s
        
			"""%(self.box1,self.box2)

		return s

def Histograms(hman,box1,box2):

	hman.h3("XYZBox1", "XYZBox1", 
		10, box1.xmin, box1.xmax,
		10, box1.ymin, box1.ymax,
		10, box1.zmin, box1.zmax)
	hman.fetch("XYZBox1").GetXaxis().SetTitle(
		"XYZ interaction box 1 (mm)")

	hman.h2("XYBox1", "XYBox1", 
		25, box1.xmin, box1.xmax,
		25, box1.ymin, box1.ymax)
	hman.fetch("XYBox1").GetXaxis().SetTitle(
		"XY interaction box 1 (mm)")

	hman.h1("XBox1", "XBox1", 
		100, box1.xmin, box1.xmax)
	hman.fetch("XBox1").GetXaxis().SetTitle(
		"X interaction box 1 (mm)")

	hman.h1("YBox1", "YBox1", 
		100, box1.ymin, box1.ymax)
	hman.fetch("YBox1").GetXaxis().SetTitle(
		"Y interaction box 1 (mm)")

	hman.h1("ZBox1", "ZBox1", 
		100, abs(box1.zmin), abs(box1.zmax))
	hman.fetch("ZBox1").GetXaxis().SetTitle(
		"Z interaction box 1 (mm)")

	hman.h3("XYZBox2", "XYZBox2", 
		10, box2.xmin, box2.xmax,
		10, box2.ymin, box2.ymax,
		10, box2.zmin, box2.zmax)
	hman.fetch("XYZBox2").GetXaxis().SetTitle(
		"XYZ interaction box 2 (mm)")

	hman.h2("XYBox2", "XYBox2", 
		25, box2.xmin, box2.xmax,
		25, box2.ymin, box2.ymax)
	hman.fetch("XYBox2").GetXaxis().SetTitle(
		"XY interaction box 2 (mm)")

	hman.h1("XBox2", "XBox2", 
		100, box2.xmin, box2.xmax)
	hman.fetch("XBox2").GetXaxis().SetTitle(
		"X interaction box 1 (mm)")

	hman.h1("YBox2", "YBox2", 
		100, box2.ymin, box2.ymax)
	hman.fetch("YBox2").GetXaxis().SetTitle(
		"Y interaction box 1 (mm)")

	hman.h1("ZBox2", "ZBox2", 
		100, abs(box2.zmin), abs(box2.zmax))
	hman.fetch("ZBox2").GetXaxis().SetTitle(
		"Z interaction box 1 (mm)")



if __name__ == '__main__':
	m = Messenger(0)
	hman =HistoManager() 
	tman =TreeManager() 

	box1 = Box(BOX1)
	box2 = Box(BOX2)
	Histograms(hman,box1,box2)

	px1 = array.array('f',[0.])
	py1 = array.array('f',[0.])
	pz1 = array.array('f',[0.])
	px2 = array.array('f',[0.])
	py2 = array.array('f',[0.])
	pz2 = array.array('f',[0.])
	

	tman.book('tpg',"photon generator tree")
	tman.addBranch('tpg','px1',px1,dim=1)
	tman.addBranch('tpg','py1',py1,dim=1)
	tman.addBranch('tpg','pz1',pz1,dim=1)
	tman.addBranch('tpg','px2',px2,dim=1)
	tman.addBranch('tpg','py2',py2,dim=1)
	tman.addBranch('tpg','pz2',pz2,dim=1)
	
		
	pg = PhotonGenerator(box1,box2,0)
	print pg

	nfail = 0
	nOK = 0
	nevents = 1000000
	nprint = 1000
	for event in range(0,nevents):
		if event%nprint == 0:
			print 'event ', event
		evt =pg.GenerateEvent()
		if evt == False:
			nfail+=1
			continue
		nOK+=1
		x1,y1,z1 = pg.X1
		x2,y2,z2 = pg.X2

		m.log(1, ' event =%d x1 =%7.2f mm, y1 =%7.2f mm, z1 =%7.2f mm '%(
					event,x1/mm,y1/mm,z1/mm))

		m.log(1, 'event =%d x2 =%7.2f mm, y2 =%7.2f mm, z2 =%7.2f mm '%(
					event,x2/mm,y2/mm,z2/mm))

		hman.fill("XYZBox1", x1/mm,y1/mm,z1/mm)
		hman.fill("XYBox1", x1/mm,y1/mm)
		hman.fill("XBox1", x1/mm)
		hman.fill("YBox1", y1/mm)
		hman.fill("ZBox1", abs(z1)/mm)

		hman.fill("XYZBox2", x2/mm,y2/mm,z2/mm)
		hman.fill("XYBox2", x2/mm,y2/mm)
		hman.fill("XBox2", x2/mm)
		hman.fill("YBox2", y2/mm)
		hman.fill("ZBox2", abs(z2)/mm)

		px1[0]=x1
		py1[0]=y1
		pz1[0]=z1

		px2[0]=x2
		py2[0]=y2
		pz2[0]=z2

		tman.fill('tpg')

	m.log(0,"nevents =%d, nfail = %d nOK =%d"%(nevents,nfail,nOK))
	pathFile = '/Users/jjgomezcadenas/Development/PETALO/WORK/histo/'
	fileName = 'GeneratePhotons.root'
	hfile = pathFile+fileName
	pathFile = '/Users/jjgomezcadenas/Development/PETALO/WORK/tree/'
	fileName = 'GeneratePhotonsTree.root'
	tfile = pathFile+fileName
	
	hman.save(file_name=hfile)
	tman.save(file_name=tfile)


	
	

