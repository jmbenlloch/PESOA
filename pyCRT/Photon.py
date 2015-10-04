import random as rnd
from BLogging import *
from LXSC import *
from ROOT import *


log =logging.getLogger("Photon")
log.setLevel(logging.INFO)
log.addHandler(ch)
deb = Debug(name="Photon",switch=0)

class SPhoton:
	def __init__(self,t0,x,y,z,tx,ty,tz,tau,lamda,n,sigmaz=1*mm):
		"""
		A scintillation photon.
		t0 = creation time stamp (time after initial interaction)
		x,y,z:  position
		tx,ty,tz: director cosines
		tau: lifetime of scintillation process that created the photon
		lamda: wavelength
		n: refraction index to wavelength 
		sigmaz: resolution in z coordinate (relevant for CTR)
		"""
		self.lamda = lamda  #photon wavelenght
		self.n = n #refraction index 
		self.sigmaz = sigmaz
		self.t0 = t0  #creation time wrt interaction time
		self.time = self.t0  #phton timestamp initially set to t0
		self.ctime = self.t0 #corrected timestamp
		self.path = 0  #accumulated path
		self.nb = 0  #number of bounces
		self.x = x
		self.y = y
		self.z = z
		self.tx = tx
		self.ty = ty
		self.tz = tz
		self.tau = tau 
		self.zr = self.z + rnd.gauss(0, self.sigmaz) # reconstructed z
		
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

	def DOI(self):
		"""
		reconstructed z coordinate (depth of interaction)
		"""
		return self.zr

	def XYZ(self):
		"""
		xyz coordinate
		"""
		return self.x,self.y,self.z

	def TX(self):
		"""
		tx cosinus tx = x/d
		"""
		return self.tx

	def TY(self):
		"""
		ty cosinus ty = y/d
		"""
		return self.ty

	def TZ(self):
		"""
		tz cosinus tz = z/d
		"""
		return self.tz

	def TXYZ(self):
		"""
		(tx,ty,tz) 
		"""
		return self.tx,self.ty,self.tz

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

	def CTime(self):
		"""
		accumulated time with respect to interaction time 
		corrected by DOI
		"""
		return self.ctime

	def N(self):
		"""
		refraction index
		"""
		return self.n

	def __str__(self):

		s= """
	      Photon:
	      x = %7.2f mm y = %7.2f mm z = %7.2f mm 
	      tx = %7.2f ty = %7.2f tz = %7.2f
	      tau = %7.4g ns,t0 = %7.4g ns 
	      time = %7.4g ns ctime = %7.4g ns path = %7.4g mm,  
	      doi =%7.4g mm
	      Lambda = %7.2f n = %7.2f
        
			"""%(self.X()/mm, self.Y()/mm, self.Z()/mm,
				self.TX(), self.TY(), self.TZ(), self.Tau()/ns,self.T0()/ns, 
				self.Time()/ns, self.CTime()/ns, self.Path()/mm,
				self.DOI()/mm, self.Lambda()/nm, self.N())

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
		"""
		test = rnd.uniform(0.,1.)

		if test > self.lxsc.PLXSC().SensorEfficiency():
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
			photon.tx = -photon.tx
			photon.ty = rnd.uniform(-1.,1.)
			photon.tz = rnd.uniform(-1.,1.)
		elif jd == 1: 
			photon.tx = rnd.uniform(-1.,1.)
			photon.ty = -photon.ty 
			photon.tz = rnd.uniform(-1.,1.)
		elif jd == 2:
			photon.tx = rnd.uniform(-1.,1.)
			photon.ty = rnd.uniform(-1.,1.)
			photon.tz = -photon.tz
		else:
			print "error: face index =",jd
			sys.exit()  


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
		x0 = photon.X()
		y0 = photon.Y()
		z0 = photon.Z()
		tx = photon.TX()
		ty = photon.TY()
		tz = photon.TZ()
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
		xx = -9999
		yy = -9999
		zz = -9999

		# if x0 == 0 or lx-x0 == 0:  #face x
		# 	return 0
		# if y0 == 0 or ly-y0 == 0:  #face y
		# 	return 1
		# if z0 == 0 or lz-z0 == 0:  #face z
		# 	return 2

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

		if -y0/ty >0 :
			D.append(-y0/ty)
			yy = 0
		elif (ly-y0)/ty >0 :
			D.append((ly-y0)/ty)
			yy = ly
		else:
			log.error(' negative distance to plane y: -y0/ty (mm) = %7.2f (ly-y0)/ty (mm) = %7.2f', 
				-(y0/ty)/mm,((ly-y0)/ty)/mm)
			sys.exit()

		if -z0/tz >0 :
			D.append(-z0/tz)
			zz = 0
			doi = photon.DOI()
		elif (lz-z0)/tz >0 :
			D.append((lz-z0)/tz)
			zz = lz
			doi = lz - photon.DOI()
		else:
			log.error(' negative distance to plane z: -z0/tz (mm) = %7.2f (lz-z0)/tz (mm) = %7.2f', 
				-(z0/tz)/mm,((lz-z0)/tz)/mm)
			sys.exit()

		d = min(D)
		jd = D.index(d)

		log.debug('next face index (1:x, 2:y, 3:z) = %d -- distance = %7.2f mm ',
				jd,d/mm)

		if jd == 0: # intersects x =0 or x =l plane
			photon.x = xx
			photon.y = y0 + d*ty
			photon.z = z0 + d*tz
		elif jd == 1:  # intersects y =0 or y =l plane
			photon.x = x0 + d*tx
			photon.y = yy
			photon.z = z0 + d*tz
		elif jd == 2:  # intersects z =0 or z =l plane
			photon.x = x0 + d*tx
			photon.y = y0 + d*ty
			photon.z = zz

		photon.path = path + d 
		photon.time = time + d*photon.n/c_light
		photon.ctime = photon.time - photon.DOI()*photon.n/c_light
		
		
		log.debug('x = %7.2f y = %7.2f z = %7.2f (in mm)', 
					photon.x/mm,photon.y/mm,photon.z/mm)
		log.debug('d = %7.2f mm --path (mm) = %7.2f -- time (ns) = %7.2f ctime (ns) = %7.2f ', 
					d,photon.path/mm,photon.time/ns, photon.ctime/ns)

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
	
		self.FTAU1 = TF1("FTAU1","exp(-x/[0])",0,10*self.lxe.tau1/ns)
		self.FTAU1.SetParameter(0,self.lxe.tau1/ns)

		self.FTAU2 = TF1("FTAU2","exp(-x/[0])",0,10*self.lxe.tau2/ns)
		self.FTAU2.SetParameter(0,self.lxe.tau2/ns)

		self.FTAU3 = TF1("FTAU3","exp(-x/[0])",0,10*self.lxe.tau3/ns)
		self.FTAU3.SetParameter(0,self.lxe.tau3/ns)
		
		self.FTAUWLS = TF1("FTAUWLS","exp(-x/[0])",0,10*self.wls.tau/ns)
		self.FTAUWLS.SetParameter(0,self.wls.tau/ns)

		self.WBox = self.lxsc.box

	def GeneratePhotonWLS(self,x=-9999,y=-9999,z=-9999,
								tx=-9999,ty=-9999,tz=-9999):
		"""
		Generates a WLS visible photon 
		"""
		self.t0 = self.FTAUWLS.GetRandom()*ns
		self.tau = self.wls.tau

		self._generatePositionAndAngles(x,y,z,tx,ty,tz)

		sp = SPhoton(self.t0,self.x,self.y,self.z,self.tx,self.ty,self.tz,self.tau,
					 self.wls.ScintillationWavelength(),
					 self.lxe.RefractionIndexBlue(),
					 self.lxsc.PLXSC().SigmaZ())
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

		sp = SPhoton(self.t0,self.x,self.y,self.z,self.tx,self.ty,self.tz,self.tau,
					 self.lxe.ScintillationWavelength(),self.lxe.RefractionIndexUV(),
					 self.lxsc.PLXSC().SigmaZ())
		return sp

	def _generatePositionAndAngles(self,x,y,z,tx,ty,tz):
		"""
		Generates the position and angles of the photon
		"""

		if x == -9999: # generate random
			self.x = rnd.random()*self.WBox.X()
		elif x > 0 and x < self.WBox.X():
			self.x = x
		else:
			log.error('error: x outside box: x (cm) = %7.2f',x/cm)
			sys.exit()
		if y == -9999:
			self.y = rnd.random()*self.WBox.Y()
		elif y > 0 and y < self.WBox.Y():
			self.y = y
		else:
			log.error('error: x outside box: y (cm) = %7.2f',y/cm)
			sys.exit()
		if z == -9999:
			self.z = rnd.random()*self.WBox.Z()
		elif z > 0 and z < self.WBox.Z():
			self.z = z
		else:
			log.error('error: x outside box: z (cm) = %7.2f',z/cm)
			sys.exit()
		if tx == -9999:
			self.tx = rnd.uniform(-1.,1.)
		elif tx >= -1 and tx <= 1:
			self.tx = tx
		else:
			log.error('error: wrong tx = %7.2f',tx)
			sys.exit()
		if ty == -9999:
			self.ty = rnd.uniform(-1.,1.)
		elif ty >= -1 and ty <= 1:
			self.ty = ty
		else:
			log.error('error: wrong ty = %7.2f',ty)
			sys.exit()
		if tz == -9999:
			self.tz = rnd.uniform(-1.,1.)
		elif tz >= -1 and tz <= 1:
			self.tz = tz
		else:
			log.error('error: wrong tz = %7.2f',tz)
			sys.exit()

	def __str__(self):
        
		s= """
        	lxe = %s
        	wls = %s
        	lxsc = %s

		"""%(self.lxe, self.wls, self.lxsc)
		return s

def GeneratePhotons():
	c1 = TCanvas( 'c1', 'CTR', 200, 10, 600, 800 )
	NPHOTONS = 10000

	lxe = LXe()
	wbox = Box(50*mm,50*mm,50*mm)
	tpb = WLS()
	sipm = SiPM()
	lxsc = LXSC(lxe,wbox,tpb,sipm,6.2*mm,[1,1,0,0,0,0])

	pg = PhotonGenerator(lxsc)
	

	htau1 = TH1F("htau1", "tau1 (ns)", 200, 0., 5*lxe.Lifetime(1)/ns)
	htau2 = TH1F("htau2", "tau2 (ns)", 200, 0., 5*lxe.Lifetime(2)/ns)
	htau3 = TH1F("htau3", "tau3 (ns)", 200, 0., 5*lxe.Lifetime(3)/ns)
	htauWLS = TH1F("htauWLS", "tauWLS (ns)", 200, 0., 5*tpb.Lifetime()/ns)

	for i in range(NPHOTONS):
		log.info('scintillation photon--> = %d',i)
				
		uvp = pg.GeneratePhotonUV(1,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		htau1.Fill(uvp.T0()/ns)
		uvp = pg.GeneratePhotonUV(2,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		htau2.Fill(uvp.T0()/ns)
		uvp = pg.GeneratePhotonUV(3,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		htau3.Fill(uvp.T0()/ns)
		bp = pg.GeneratePhotonWLS(x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		htauWLS.Fill(bp.T0()/ns)

	c1.Divide(2,2)
	c1.cd(1)
	gPad.SetLogy()
	htau1.Draw()
	
	c1.cd(2)
	gPad.SetLogy()
	htau2.Draw()

	c1.cd(3)
	gPad.SetLogy()
	htau3.Draw()

	c1.cd(4)
	gPad.SetLogy()
	htauWLS.Draw()

	c1.Show()
	wait()


def ScintillationEvent(tauIndex):
	"""
	Generate and propagate photons to the readout faces
	tauIndex =1,2,3 for one of the three lifetimes of Xe  
	"""

	c1 = TCanvas( 'c1', 'CTR', 200, 10, 600, 800 )
	NPHOTONS = 20000
	MAX_BOUNCES=100

	lxe = LXe()
	# plxsc = PLXSC(sigmax=1e-3*mm, sigmay=1e-3*mm, sigmaz=1e-3*mm, sigma0=0.05,
 #               sensorEff=1.0, wlsEff=1.0, uvRef=1.0, wlsRef=1.0)

	plxsc = PLXSC(sigmax=1*mm, sigmay=1*mm, sigmaz=1*mm, sigma0=0.05,
               sensorEff=0.9, wlsEff=0.9, uvRef=0.95, wlsRef=0.98)

	wbox = Box(50*mm,50*mm,50*mm)
	tpb = WLS()
	sipm = SiPM()
	imask = [1,0,0]
	lxsc = LXSC(lxe,wbox,tpb,plxsc,sipm,6.2*mm,imask)
	
	pg = PhotonGenerator(lxsc)
	pt = PhotonTransport(lxsc)
	
	htau = TH1F("htau", "tau (ns)", 200, 0., 5*lxe.Lifetime(tauIndex)/ns)
	htime = TH1F("htime", "photon time (ns)", 200, 0., 5*lxe.Lifetime(tauIndex)/ns)
	hctime = TH1F("hctime", "photon time (ns)", 200, 0., 5*lxe.Lifetime(tauIndex)/ns)
	hpath = TH1F("hpath", "photon path (mm)", 100, 0., 10*wbox.X())
	hb = TH1F("hb", "number of bounces", 50, 0., 50.)
	
	for i in range(NPHOTONS):
		log.info('scintillation photon--> = %d',i)
		log.debug('tau = %7.2f ns x = y = z = %7.2f mm',lxe.Lifetime(tauIndex)/ns,(wbox.X()/2)/mm)

		#Generate VUV photon
		uvp = pg.GeneratePhotonUV(tauIndex,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		htau.Fill(uvp.T0()/ns)

		log.debug('Generating UV photon = %s',uvp)
		deb.Wait()

		
		n=0
		while n < MAX_BOUNCES: #propagate up to the max number of bounces
			fi = pt.Step(uvp)  #step the photon to the next face

			log.debug('face	index  = %d',fi)

			
			if pt.TestInstrumentedFace(fi) == 1: #hitting instrumented face, photon is absorbed
				log.debug('photon hits instrumented face index  = %d',fi)

				#Efficiency of detection by sensor plane

				if pt.TestSensorEfficiency() == True:
					log.debug('photon is detected by sensors')
					
					#histogram path and time
					htime.Fill(uvp.Time()/ns)
					hpath.Fill(uvp.Path()/mm)

					#Sensor Response goes here.
					hctime.Fill(uvp.CTime()/ns)
					break

				else:
					log.debug('photon not detected')
					break
				

			else: #is VUV photon absorbed by Teflon?
				test = pt.TestReflectivityUV()
					
				if test == False:
					log.debug('photon absorbed')
					break
				else: # Lambertian reflection
					log.debug('Lambertian reflection in face index = %d',fi)
					pt.Lambert(uvp,fi) 
					hb.Fill(n)
			n+=1
			
	
	c1.Divide(2,2)
	c1.cd(1)
	gPad.SetLogy()
	htau.Draw()
	
	c1.cd(2)
	gPad.SetLogy()
	htime.Draw()

	c1.cd(3)
	gPad.SetLogy()
	hctime.Draw()

	c1.cd(4)
	gPad.SetLogy()
	hpath.Draw()

	c1.Show()
	wait()


if __name__ == '__main__':
	#GeneratePhotons()
	
	ScintillationEvent(1)


