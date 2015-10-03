import random as rnd
import logging 
from LXSC import *
from ROOT import *

ch = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)

log =logging.getLogger("Photon")
log.setLevel(logging.DEBUG)
log.addHandler(ch)

class SPhoton:
	def __init__(self,t0,x,y,z,tx,ty,tz,tau):
		"""
		A scintillation photon.
		"""
		self.lamda = 172*nm  #photon is created as VUV photon
		self.n = 1.55 #refraction index for VUV
		self.t0 = t0  #creation time wrt interaction time
		self.time = self.t0  #phton timestamp initially set to t0
		self.path = 0  #accumulated path
		self.nb = 0  #number of bounces
		self.x = x
		self.y = y
		self.z = z

		self.tx = tx
		self.ty = ty
		self.tz = tz
		self.tau = tau 
		
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

	def __str__(self):

		s= """
	      Photon:
	      x = %7.2f mm y = %7.2f mm z = %7.2f mm
	      tx = %7.2f ty = %7.2f tz = %7.2f
	      tau = %7.4g ns,t0 = %7.4g ns 
	      time = %7.4g ns path = %7.4g mm, 
	      Lambda = %7.2f
        
			"""%(self.X()/mm, self.Y()/mm, self.Z()/mm,
				self.TX(), self.TY(), self.TZ(), self.Tau()/ns,self.T0()/ns, 
				self.Time()/ns, self.Path()/mm, self.Lambda()/nm )

		return s

class PhotonTransport:
	"""
	Transports the photon to the next face in the box.
	The box is defined by six planes: 
	plane 0 x =0, (x = lx, where lx is the length in x of the box )
	plane 1 y =0, (y = ly, where ly is the length in y of the box )
	plane 2 z =0, (z = lz, where lz is the length in z of the box )
	"""
	def __init__(self,box):
		self.box = box

	def  lambert(self,photon,jd):
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


	def step(self,photon):
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

		lx = self.box.X()
		ly = self.box.Y()
		lz = self.box.Z()
		
		log.debug('x0 = %7.2f y0 = %7.2f z0 = %7.2f', x0/mm,y0/mm,z0/mm)
		log.debug('tx = %7.2f ty = %7.2f tz = %7.2f', tx,ty,tz)

		D=[]
		xx = -9999
		yy = -9999
		zz = -9999

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
		elif (lz-z0)/tz >0 :
			D.append((lz-z0)/tz)
			zz = ly
		else:
			log.error(' negative distance to plane z: -z0/tz (mm) = %7.2f (lz-z0)/tz (mm) = %7.2f', 
				-(z0/tz)/mm,((lz-z0)/tz)/mm)
			sys.exit()

		# print D
		# print min(D)
		# print D.index(min(D))

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

		path = photon.path + d
		photon.path = path 
		time = photon.time + photon.path/C
		photon.time = time

		log.debug('x = %7.2f y = %7.2f z = %7.2f (in mm)', 
					photon.x/mm,photon.y/mm,photon.z/mm)
		log.debug('path (mm) = %7.2f time (ns) = %7.2f ', 
					photon.path/mm,photon.time/ns)

		return jd


class PhotonGenerator:
	"""
	Generates a Xenon scintillation photon
	"""
	def __init__(self,WBox,tau1 = 2.2*ns,tau2 = 27*ns,tau3 = 45*ns):
		self.tau1 = tau1
		self.tau2 = tau2
		self.tau3 = tau3
	
		self.FTAU1 = TF1("FTAU1","exp(-x/[0])",0,100)
		self.FTAU1.SetParameter(0,self.tau1/ns)

		self.FTAU2 = TF1("FTAU2","exp(-x/[0])",0,300)
		self.FTAU2.SetParameter(0,self.tau2/ns)

		self.FTAU3 = TF1("FTAU3","exp(-x/[0])",0,500)
		self.FTAU3.SetParameter(0,self.tau3/ns)

		self.WBox = WBox

	def GeneratePhoton(self,tau=2.2*ns,x=-9999,y=-9999,z=-9999,
										tx=-9999,ty=-9999,tz=-9999):
		"""
		Generates tha actual photon 
		"""
		if tau == self.tau1:
			self.t0 = self.FTAU1.GetRandom()*ns
		elif tau == self.tau2:
			self.t0 = self.FTAU2.GetRandom()*ns
		elif tau == self.tau3:
			self.t0 = self.FTAU3.GetRandom()*ns
		else:
			log.error(' wrong lifetime: it must be %7.2f ns or %7.2f ns or %7.2f ns ',
				self.tau1/ns,self.tau2/ns,self.tau3/ns)
			sys.exit()
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

		sp = SPhoton(self.t0,self.x,self.y,self.z,self.tx,self.ty,self.tz,tau)
		return sp

if __name__ == '__main__':

	c1 = TCanvas( 'c1', 'CTR', 200, 10, 600, 800 )
	
	NPHOTONS = 1000
	TAU1 = 2.2*ns
	TAU2 = 27*ns
	TAU3 = 45*ns
	LX = 5*cm
	LY = 5*cm
	LZ = 5*cm

	
	htau1 = TH1F("htau1", "tau1 (ns)", 200, 0., 5*TAU1/ns)
	htau2 = TH1F("htau2", "tau2 (ns)", 200, 0., 5*TAU2/ns)
	htau3 = TH1F("htau3", "tau3 (ns)", 200, 0., 5*TAU3/ns)
	

	wbox = Box(LX,LY,LZ)
	pt = PhotonTransport(wbox)
	pg = PhotonGenerator(wbox)

	for i in range(NPHOTONS):
		log.info('scintillation photon--> = %d',i)
		
		uvp = pg.GeneratePhoton(tau=TAU1,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		htau1.Fill(uvp.T0()/ns)
		uvp = pg.GeneratePhoton(tau=TAU2,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		htau2.Fill(uvp.T0()/ns)
		uvp = pg.GeneratePhoton(tau=TAU3,x=wbox.X()/2,y=wbox.Y()/2,z=wbox.Z()/2)
		htau3.Fill(uvp.T0()/ns)


	c1.Divide(1,3)
	c1.cd(1)
	gPad.SetLogy()
	htau1.Draw()
	
	c1.cd(2)
	gPad.SetLogy()
	htau2.Draw()

	c1.cd(3)
	gPad.SetLogy()
	htau3.Draw()

	c1.Show()
	wait()
