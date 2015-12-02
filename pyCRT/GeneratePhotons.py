import random as rnd
from BLogging import *
from LXSC import *
from THisto import *

log =logging.getLogger("PhotonGenerator")
log.setLevel(logging.INFO)
log.addHandler(ch)
deb = Switch(name="PhotonGenerator",switch=0)

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
	def __init__(self,boxCoord1,boxCoord2,nevents):
		"""
		BOX = the boxes
		nevents = number of events to generate
		 
		
		"""
		self.box1 = Box(boxCoord1)
		self.box2 = Box(boxCoord2)  
		self.nevents = nevents
		
		
	def GenerateMomentum(self):
		"""
		generate (px,py,pz) for photon 1
		"""
		return Point3D(self.xyzt.x,self.xyzt.y,self.xyzt.z)

	def Time(self):
		"""
		actual time with respect to interaction time 
		"""
		return self.xyzt.t

	def DirectionCosines(self):
		"""
		actual direction cosines
		"""
		return self.txyz

	def T0(self):
		"""
		Creation time with respect to interaction time 
		"""
		return self.xyzt0.t

	def Position0(self):
		"""
		Position  with respect to interaction time 
		"""
		return Point3D(self.xyzt0.x,self.xyzt0.y,self.xyzt0.z)

	def NBounces(self):
		"""
		number of bounces
		"""
		return self.nb

	def Lambda(self):
		"""
		wavelenght: (172 or 420)
		"""
		return self.lamda  

	def Tau(self):
		"""
  		Decay lifetime of the process that generated this photon
		"""
		return self.tau

	def Path(self):
		"""
		accumulated path  
		"""
		return self.path
	
	def N(self):
		"""
		refraction index
		"""
		return self.n

	def __str__(self):

		s= """
		  initial coordinates (x,y,z,t) = %s
		  actual coordinates (x,y,z,t) = %s
		  actual cosines (tx,ty,tz) = %s
	      tau = %7.4g ns, 
	      path = %7.4g mm,  
	      Lambda = %7.2f n = %7.2f
        
			"""%(self.xyzt0,self.xyzt, self.txyz,
				self.Tau()/ns, 
				self.Path()/mm,
				self.Lambda()/nm, self.N())

		return s

class Photoelectron:
	"""
	Represents a photoelectron:
	Takes the initial photon instance,  and the sipm number.
	Only the photon is set at instance time 
	"""
	def __init__(self,photon):
		self.photon = photon
		self.ctime = 0
		self.crtime = 0
		self.nsipm = 0

	def __str__(self):

		s= """
		  Photon = %s
		  ctime = %7.2f ns crtime = %7.2f ns 
		  SiPM number = %d
        
			"""%(self.photon, self.ctime, self.crtime, self.nsipm)

		return s


class PhotonTransport:
	"""
	Transports the photon to the next face in the LXSC box.
	The box is defined by six planes: 
	plane 0 x =0, (x = lx, where lx is the length in x of the box )
	plane 1 y =0, (y = ly, where ly is the length in y of the box )
	plane 2 z =0, (z = lz, where lz is the length in z of the box )
	"""
	def __init__(self,lxsc,detEff=0.9,lxeReflectivity=0.95):
		self.lxsc = lxsc
		self.box = self.lxsc.Box()
		self.imask = self.lxsc.instMask
		

	def TestInstrumentedFace(self,fi):
		"""
		Tests whether the photon hits an instrumented face as specified
		by the mask imask.
		fi = 0 for x =0 or x=l; 1 for y=0 or y = l; 2 for z=0 or z= l
		a mask [1,0,0] means that x is instrumented
		a mask [1,1,1] means that x,y,z are instrumented
		thus if instrumented imask[fi]=1
		"""

		return self.imask[fi]

	def TestSensorEfficiency(self):
		"""
		Tests whether photon is detected by sensors
		(geomeatrical coverage)
		"""
		test = rnd.uniform(0.,1.)

		if test > self.lxsc.SensorCoverage():
			return False
		else:
			return True

	def TestSensorPDE(self):
		"""
		Tests PDE
		"""
		test = rnd.uniform(0.,1.)

		if test > self.lxsc.Sensor().PDE():
			return False
		else:
			return True


	def TestReflectivityUV(self):
		"""
		Tests whether UV photon reflects in the box wall or is absorbed 
		"""
		test = rnd.uniform(0.,1.)

		if test > self.lxsc.PLXSC().ReflectivityUV():
			return False
		else:
			return True 

	def TestWLSEfficiency(self):
		"""
		Tests whether UV photon emits a WLS photon  
		"""
		test = rnd.uniform(0.,1.)

		if test > self.lxsc.PLXSC().WLSEfficiency():
			return False
		else:
			return True

	def TestReflectivityWLS(self):
		"""
		Tests whether wls (blue) photon reflects in the box wall or is absorbed 
		"""
		test = rnd.uniform(0.,1.)

		if test > self.lxsc.PLXSC().ReflectivityWLS():
			return False
		else:
			return True 

					
	def Lambert(self,photon,jd):
		"""
		Lambertian reflection
		changes the direction of the incoming cosinus
		random generation of the other two cosinus. 
		"""
		if jd == 0:
			photon.txyz.x = -photon.txyz.x
			photon.txyz.y = rnd.uniform(-1.,1.)
			photon.txyz.z = rnd.uniform(-1.,1.)
		elif jd == 1: 
			photon.txyz.x = rnd.uniform(-1.,1.)
			photon.txyz.y = -photon.txyz.y
			photon.txyz.z = rnd.uniform(-1.,1.)

		elif jd == 2:
			photon.txyz.x = rnd.uniform(-1.,1.)
			photon.txyz.y = rnd.uniform(-1.,1.)
			photon.txyz.z = -photon.txyz.z 
		else:
			print "error: face index =",jd
			sys.exit()  

	def ComputePath(self,D,lx,x0,tx):
		"""
		to step the photon to the plane x=0 (x=lx) one checks the equations:
		x = 0 = x0 + d*tx --> d = -x0/tx
		x = lx = x0 + d*tx --> d = (lx-x0)/tx

		"""
		xx = -9999

		if -x0/tx >0 :
			D.append(-x0/tx)
			xx = 0
		elif (lx-x0)/tx >0 :
			D.append((lx-x0)/tx)
			xx = lx
		else:
			log.error(' negative distance to plane x: -x0/tx (mm) = %7.2f (lx-x0)/tx (mm) = %7.2f', 
				-(x0/tx)/mm,((lx-x0)/tx)/mm)
			sys.exit()
		return xx

	def Step(self,photon):
		"""
		Steps to the closer face of the box.
		the photon is defined by coordinates:
		x = x0 + d*tx
		y = y0 + d*ty
		z = z0 + d*tz

		to step the photon to the plane x=0 (x=lx) one checks the equations:
		x = 0 = x0 + d*tx --> d = -x0/tx
		x = lx = x0 + d*tx --> d = (lx-x0)/tx

		and takes the distance that is positive. Repeat for the other two 
		coordinates and take the shortest distance.

		"""
		x0 = photon.Position().X()
		y0 = photon.Position().Y()
		z0 = photon.Position().Z()
		tx = photon.DirectionCosines().X()
		ty = photon.DirectionCosines().Y()
		tz = photon.DirectionCosines().Z()
		path = photon.Path()
		time = photon.Time()

		lx = self.box.X()
		ly = self.box.Y()
		lz = self.box.Z()

		log.debug('+++step+++')
		log.debug('x0 = %7.2f mm y0 = %7.2f mm z0 = %7.2f mm', x0/mm,y0/mm,z0/mm)
		log.debug('tx = %7.2f ty = %7.2f tz = %7.2f', tx,ty,tz)
		log.debug('path = %7.2f mm time = %7.2f ns ', path,time)
		
		deb.Wait()

		D=[]
		
		xx = self.ComputePath(D,lx,x0,tx)
		yy = self.ComputePath(D,ly,y0,ty)
		zz = self.ComputePath(D,lz,z0,tz)

		d = min(D)
		jd = D.index(d)

		log.debug('next face index (1:x, 2:y, 3:z) = %d -- distance = %7.2f mm ',
				jd,d/mm)

		if jd == 0: # intersects x =0 or x =l plane
			photon.xyzt.x = xx
			photon.xyzt.y = y0 + d*ty
			photon.xyzt.z = z0 + d*tz
		elif jd == 1:  # intersects y =0 or y =l plane
			photon.xyzt.x = x0 + d*tx
			photon.xyzt.y = yy
			photon.xyzt.z = z0 + d*tz
		elif jd == 2:  # intersects z =0 or z =l plane
			photon.xyzt.x = x0 + d*tx
			photon.xyzt.y = y0 + d*ty
			photon.xyzt.z = zz

		photon.path = path + d 
		photon.xyzt.t = time + d*photon.N()/c_light
		
		log.debug('step: Photon at xyzt = %s',photon.xyzt)
		log.debug('step: Photon at xyzt0 = %s',photon.xyzt0)
		log.debug('d = %7.2f mm --path (mm) = %7.2f -- time (ps) = %7.2f ', 
					d,photon.Path()/mm,photon.Time()/ps)

		deb.Wait()

		return jd


class PhotonGenerator:
	"""
	Generates a Xenon scintillation photon
	"""
	def __init__(self,lxsc):
		self.lxsc = lxsc
		self.lxe = self.lxsc.LXe()
		self.wls = self.lxsc.WLS()
		self.WBox = self.lxsc.box
		self.position = Point3D()
		self.cosines = Point3D()
	
		self.FTAU1 = TF1("FTAU1","exp(-x/[0])",0,10*self.lxe.tau1/ns)
		self.FTAU1.SetParameter(0,self.lxe.tau1/ns)

		self.FTAU2 = TF1("FTAU2","exp(-x/[0])",0,10*self.lxe.tau2/ns)
		self.FTAU2.SetParameter(0,self.lxe.tau2/ns)

		self.FTAU3 = TF1("FTAU3","exp(-x/[0])",0,10*self.lxe.tau3/ns)
		self.FTAU3.SetParameter(0,self.lxe.tau3/ns)
		
		self.FTAUWLS = TF1("FTAUWLS","exp(-x/[0])",0,10*self.wls.tau/ns)
		self.FTAUWLS.SetParameter(0,self.wls.tau/ns)

		

	def GeneratePhotonWLS(self,x=-9999,y=-9999,z=-9999,
								tx=-9999,ty=-9999,tz=-9999):
		"""
		Generates a WLS visible photon 
		"""
		t = self.FTAUWLS.GetRandom()*ns
		
		self._generatePositionAndAngles(x,y,z,tx,ty,tz)
		
		sp = SPhoton(t,self.position,self.cosines,self.wls.tau,
					 self.wls.ScintillationWavelength(),
					 self.lxe.RefractionIndexBlue())
					 
		return sp

	def GeneratePhotonUV(self,i,x=-9999,y=-9999,z=-9999,
								tx=-9999,ty=-9999,tz=-9999):
		"""
		Generates a VUV photon 
		"""
		t = 0
		tau = 0
		if i==1:
			t = self.FTAU1.GetRandom()*ns
			tau = self.lxe.tau1
		elif i == 2:
			t = self.FTAU2.GetRandom()*ns
			tau = self.lxe.tau2
		elif i == 3:
			t = self.FTAU3.GetRandom()*ns
			tau = self.lxe.tau3
		else:
			log.error(' lifetime index must be 1,2,3, index = %d ',i)
			sys.exit()

		log.debug('In GeneratePhotonUV: t0 = %7.2f ns ',
			t/ns)
		log.debug('In GeneratePhotonUV: x = %7.2f mm y = %7.2f mm z = %7.2f mm',
			x/mm,y/mm,z/mm)
		log.debug('In GeneratePhotonUV: tx = %7.2f mm ty = %7.2f mm tz = %7.2f mm',
			tx/mm,ty/mm,tz/mm)
		self._generatePositionAndAngles(x,y,z,tx,ty,tz)

		log.debug('In GeneratePhotonUV: position = %s',self.position)
		log.debug('In GeneratePhotonUV: cosines = %s',self.cosines)

			
		sp = SPhoton(t,self.position,self.cosines,tau,
					 self.lxe.ScintillationWavelength(),self.lxe.RefractionIndexUV())

		log.debug("generated photon =%s",sp)
		deb.Wait()
					 
		return sp

	def _generatePoint(self,x):
		"""
		Generates one point 
		"""
		xx=0

		if x == -9999: # generate random
			xx= rnd.random()*self.WBox.X()
		elif x > 0 and x < self.WBox.X():
			xx = x
		else:
			log.error('error: x outside box: x (cm) = %7.2f',x/cm)
			sys.exit()
		return xx

	def _generateCosine(self,x):
		"""
		Generates one point 
		"""
		xx=0

		if x == -9999: # generate random
			xx= rnd.uniform(-1.,1.)
		elif x >= -1 and x <= 1:
			xx = x
		else:
			log.error('error: tx outside box: tx  = %7.2f',x)
			sys.exit()
		return xx

	def _generatePositionAndAngles(self,x,y,z,tx,ty,tz):
		"""
		Generates the position and angles of the photon
		"""
		self.position.x = self._generatePoint(x)
		self.position.y = self._generatePoint(y)
		self.position.z = self._generatePoint(z)
		self.cosines.x = self._generateCosine(tx)
		self.cosines.y = self._generateCosine(ty)
		self.cosines.z = self._generateCosine(tz)

		
	def __str__(self):
        
		s= """
        	lxe = %s
        	wls = %s
        	lxsc = %s

		"""%(self.lxe, self.wls, self.lxsc)
		return s

def Histos(th):
	lxe = LXe()
	tpb = WLS()
	th.BookH1("htau1", "tau1 (ns)", 200, 0., 5*lxe.Lifetime(1)/ns)
	th.BookH1("htau2", "tau2 (ns)", 200, 0., 5*lxe.Lifetime(2)/ns)
	th.BookH1("htau3", "tau3 (ns)", 200, 0., 5*lxe.Lifetime(3)/ns)
	th.BookH1("htauWLS", "tauWLS (ns)", 200, 0., 5*tpb.Lifetime()/ns)

def GeneratePhotons():
	#c1 = TCanvas( 'c1', 'CTR', 200, 10, 600, 800 )
	th = THisto()
	Histos(th)
    
	NPHOTONS = 10000

	lxe = LXe()
	wbox = Box(50*mm,50*mm,50*mm)
	tpb = WLS(name="TPB", lamda=420*nm, tau=1*ns)
	sipm = Sensor()
	lxsc = LXSC(lxe,wbox,tpb,sipm,6.2*mm,[0,0,1])

	pg = PhotonGenerator(lxsc)
	

	for i in range(NPHOTONS):
		log.info('scintillation photon--> = %d',i)
				
		uvp = pg.GeneratePhotonUV(1,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		th.FillH1("htau1",uvp.T0()/ns)
		
		uvp = pg.GeneratePhotonUV(2,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		
		th.FillH1("htau2",uvp.T0()/ns)
		uvp = pg.GeneratePhotonUV(3,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		#
		th.FillH1("htau3",uvp.T0()/ns)
		bp = pg.GeneratePhotonWLS(x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		th.FillH1("htauWLS",bp.T0()/ns)
		

	th.DrawList(["htau1","htau2","htau3","htauWLS"],xd=2,yd=2)




if __name__ == '__main__':
	GeneratePhotons()
	
	

