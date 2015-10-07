import random as rnd
from BLogging import *
from LXSC import *
from THisto import *

log =logging.getLogger("Photon")
log.setLevel(logging.INFO)
log.addHandler(ch)
deb = Switch(name="Photon",switch=0)

class SPhoton:
	def __init__(self,t0,xyz,txyz,tau=2.2*ns,lamda=172*nm,n=1.55):
		"""
		A scintillation photon.
		t0 = creation time stamp (time after initial interaction)
		xyz:  a Point3D describing position
		txyz: a Point 3D describing director cosines
		tau: lifetime of scintillation process that created the photon
		lamda: wavelength
		n: refraction index to wavelength 
		
		"""
		self.lamda = lamda  #photon wavelenght
		self.n = n #refraction index 
		self.tau = tau 
		self.t0 = t0  #creation time wrt interaction time
		self.time = self.t0  #phton timestamp initially set to t0
		self.path = 0  #accumulated path
		self.nb = 0  #number of bounces
		self.xyz = xyz 
		self.txyz = txyz

		
	def Position(self):
		"""
		Position
		"""
		return self.xyz

	def DirectionCosines(self):
		"""
		cosines
		"""
		return self.txyz

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

	def T0(self):
		"""
		Creation time with respect to interaction time 
		"""
		return self.t0

	def Time(self):
		"""
		accumulated time with respect to interaction time 
		"""
		return self.time

	def N(self):
		"""
		refraction index
		"""
		return self.n

	def __str__(self):

		s= """
	     
	      t0 = %7.4g ns time = %7.4g ns
	      x = %7.2f mm y = %7.2f mm z = %7.2f mm 
	      tx = %7.2f ty = %7.2f tz = %7.2f
	      tau = %7.4g ns, 
	      path = %7.4g mm,  
	      Lambda = %7.2f n = %7.2f
        
			"""%(self.T0()/ns,self.Time()/ns,
				self.Position().X()/mm, self.Position().Y()/mm, self.Position().Z()/mm,
				self.DirectionCosines().TX(), self.DirectionCosines().TY(), 
				self.DirectionCosines().TZ(), 
				self.Tau()/ns, 
				self.Path()/mm,
				self.Lambda()/nm, self.N())

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
			photon.DirectionCosines().x = -photon.DirectionCosines().X()
			photon.DirectionCosines().y = rnd.uniform(-1.,1.)
			photon.DirectionCosines().z = rnd.uniform(-1.,1.)
		elif jd == 1: 
			photon.DirectionCosines().x = rnd.uniform(-1.,1.)
			photon.DirectionCosines().y = -photon.DirectionCosines().Y()
			photon.DirectionCosines().z = rnd.uniform(-1.,1.)

		elif jd == 2:
			photon.DirectionCosines().x = rnd.uniform(-1.,1.)
			photon.DirectionCosines().y = rnd.uniform(-1.,1.)
			photon.DirectionCosines().z = -photon.DirectionCosines().Z() 
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
		
		xx = self.ComputePath(self,D,lx,x0,tx)
		yy = self.ComputePath(self,D,ly,y0,ty)
		zz = self.ComputePath(self,D,lz,z0,tz)

		d = min(D)
		jd = D.index(d)

		log.debug('next face index (1:x, 2:y, 3:z) = %d -- distance = %7.2f mm ',
				jd,d/mm)

		if jd == 0: # intersects x =0 or x =l plane
			photon.Position().x = xx
			photon.Position().y = y0 + d*ty
			photon.Position().z = z0 + d*tz
		elif jd == 1:  # intersects y =0 or y =l plane
			photon.Position().x = x0 + d*tx
			photon.Position().y = yy
			photon.Position().z = z0 + d*tz
		elif jd == 2:  # intersects z =0 or z =l plane
			photon.Position().x = x0 + d*tx
			photon.Position().y = y0 + d*ty
			photon.Position().z = zz

		photon.path = path + d 
		photon.time = time + d*photon.N()/c_light
		
		
		log.debug('step: x = %7.2f y = %7.2f z = %7.2f (in mm)', 
					photon.Position().X()/mm,photon.Position().Y()/mm,
					photon.Position().Z()/mm)
		log.debug('d = %7.2f mm --path (mm) = %7.2f -- time (ns) = %7.2f ', 
					d,photon.Path()/mm,photon.Time()/ns)

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
		self.t0 = self.FTAUWLS.GetRandom()*ns
		self.tau = self.wls.tau
		
		self._generatePositionAndAngles(x,y,z,tx,ty,tz)

		sp = SPhoton(self.t0,self.position,self.cosines,self.tau,
					 self.wls.ScintillationWavelength(),
					 self.lxe.RefractionIndexBlue())
					 
		return sp

	def GeneratePhotonUV(self,i,x=-9999,y=-9999,z=-9999,
								tx=-9999,ty=-9999,tz=-9999):
		"""
		Generates a VUV photon 
		"""
		if i==1:
			self.t0 = self.FTAU1.GetRandom()*ns
			self.tau = self.lxe.tau1
		elif i == 2:
			self.t0 = self.FTAU2.GetRandom()*ns
			self.tau = self.lxe.tau2
		elif i == 3:
			self.t0 = self.FTAU3.GetRandom()*ns
			self.tau = self.lxe.tau3
		else:
			log.error(' lifetime index must be 1,2,3, index = %d ',i)
			sys.exit()

		self._generatePositionAndAngles(x,y,z,tx,ty,tz)

		sp = SPhoton(self.t0,self.position,self.cosines,self.tau,
					 self.lxe.ScintillationWavelength(),self.lxe.RefractionIndexUV())
					 
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

	def _generatePositionAndAngles(self,x,y,z,tx,ty,tz):
		"""
		Generates the position and angles of the photon
		"""
		self.position.x = self._generatePoint(x)
		self.position.y = self._generatePoint(y)
		self.position.z = self._generatePoint(z)
		self.cosines.x = self._generatePoint(tx)
		self.cosines.y = self._generatePoint(ty)
		self.cosines.z = self._generatePoint(tz)

		
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
	

	# htau1 = TH1F("htau1", "tau1 (ns)", 200, 0., 5*lxe.Lifetime(1)/ns)
	# htau2 = TH1F("htau2", "tau2 (ns)", 200, 0., 5*lxe.Lifetime(2)/ns)
	# htau3 = TH1F("htau3", "tau3 (ns)", 200, 0., 5*lxe.Lifetime(3)/ns)
	# htauWLS = TH1F("htauWLS", "tauWLS (ns)", 200, 0., 5*tpb.Lifetime()/ns)

	for i in range(NPHOTONS):
		log.info('scintillation photon--> = %d',i)
				
		uvp = pg.GeneratePhotonUV(1,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		th.FillH1("htau1",uvp.T0()/ns)
		#htau1.Fill(uvp.T0()/ns)
		uvp = pg.GeneratePhotonUV(2,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		#htau2.Fill(uvp.T0()/ns)
		th.FillH1("htau2",uvp.T0()/ns)
		uvp = pg.GeneratePhotonUV(3,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		#htau3.Fill(uvp.T0()/ns)
		th.FillH1("htau3",uvp.T0()/ns)
		bp = pg.GeneratePhotonWLS(x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		th.FillH1("htauWLS",bp.T0()/ns)
		#htauWLS.Fill(bp.T0()/ns)

	th.DrawList(["htau1","htau2","htau3","htauWLS"],xd=2,yd=2)

	# c1.Divide(2,2)
	# c1.cd(1)
	# gPad.SetLogy()
	# htau1.Draw()
	
	# c1.cd(2)
	# gPad.SetLogy()
	# htau2.Draw()

	# c1.cd(3)
	# gPad.SetLogy()
	# htau3.Draw()

	# c1.cd(4)
	# gPad.SetLogy()
	# htauWLS.Draw()

	# c1.Show()
	# wait()




if __name__ == '__main__':
	GeneratePhotons()
	
	

