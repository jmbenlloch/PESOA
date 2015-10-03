from math import *
import random as rnd
import numpy as np
#import pylab as P
from LXSC import *
from ROOT import *

LXE_REFLECTIVITY = 0.95

NG = 30660
NG_tau1 = 3200 
NG_tau2 = 5690 
NG_tau3 = NG- NG_tau1 - NG_tau2 

TAU1 = 2.2*ns
TAU2 = 27*ns
TAU3 = 45*ns

C = 30*cm/ns

SIPM_PR = 0.8*ns
SIPM_PD = 1.3*ns # fast pulse
SIPM_PDE = 0.2
SIPM_J = 0.1

DETECTION_EFF=0.9 #overall efficiency of detection plane

TAUMAX = 3*TAU1
TBINS = 0.01*ns
NBINS = (int) (TAUMAX/TBINS)

TLEVEL = 1 # in pes

MAX_BOUNCES = 100


NPHOTONS = NG_tau1
NEVENTS=100
TAU = TAU1
LX = 5*cm
LY = 5*cm
LZ = 5*cm

DEBUG=0
	

# FTAU1 = TF1("FTAU1","exp(-x/[0])",0,50)
# FTAU1.SetParameter(0,2.2/ns)

# FPULSE = TF1("FPULSE","(exp(-x/[0])-exp(-x/[1]))/([0]-[1])",0,10)
# FPULSE.SetParameter(0,SIPM_PD/ns)
# FPULSE.SetParameter(1,SIPM_PR/ns)

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
	Transports the photon to the next face in the box
	"""
	def __init__(self,box):
		self.box = box
		
	def  lambert(self,photon,jd):
		"""
		Lambertian reflection
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
		Steps to the closer face
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

		# print "lx,ly,lz in mm",lx/mm,ly/mm,lz/mm
		# print "x0,y0,z0 in mm",x0/mm,y0/mm,z0/mm
		# print "tx,ty,tz ",tx,ty,tz

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
			print "error: -x0/tx =",-x0/tx," (lx-x0)/tx =",(lx-x0)/tx
			sys.exit()

		if -y0/ty >0 :
			D.append(-y0/ty)
			yy = 0
		elif (ly-y0)/ty >0 :
			D.append((ly-y0)/ty)
			yy = ly
		else:
			print "error: -y0/ty =",-y0/ty," (ly-y0)/ty =",(ly-y0)/ty
			sys.exit()

		if -z0/tz >0 :
			D.append(-z0/tz)
			zz = 0
		elif (lz-z0)/tz >0 :
			D.append((lz-z0)/tz)
			zz = ly
		else:
			print "error: -z0/tz =",-z0/tz," (lz-z0)/tz =",(lz-z0)/tz
			sys.exit()

		# print D
		# print min(D)
		# print D.index(min(D))

		d = min(D)
		jd = D.index(d)

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
		#print "x,y,z in mm",photon.x/mm,photon.y/mm,photon.z/mm
		#print "path in mm, time in ns",photon.path/mm,photon.time/mm

		return jd


class PhotonGenerator:
	"""
	Generates a photon
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
			print "wrong lifetime chosen"
			sys.exit()
		if x == -9999:
			self.x = rnd.random()*self.WBox.X()
		elif x > 0 and x < self.WBox.X():
			self.x = x
		else:
			print "error: x outside box: x (cm) = ",x/cm
			sys.exit()
		if y == -9999:
			self.y = rnd.random()*self.WBox.Y()
		elif y > 0 and y < self.WBox.Y():
			self.y = y
		else:
			print "error: y outside box: y (cm) = ",y/cm
			sys.exit()
		if z == -9999:
			self.z = rnd.random()*self.WBox.Z()
		elif z > 0 and z < self.WBox.Z():
			self.z = z
		else:
			print "error: z outside box: z (cm) = ",z/cm
			sys.exit()
		if tx == -9999:
			self.tx = rnd.uniform(-1.,1.)
		elif tx >= -1 and tx <= 1:
			self.tx = tx
		else:
			print "error: wrong tx  = ",tx
			sys.exit()
		if ty == -9999:
			self.ty = rnd.uniform(-1.,1.)
		elif ty >= -1 and ty <= 1:
			self.ty = ty
		else:
			print "error: wrong ty  = ",ty
			sys.exit()
		if tz == -9999:
			self.tz = rnd.uniform(-1.,1.)
		elif tz >= -1 and tz <= 1:
			self.tz = tz
		else:
			print "error: wrong tz  = ",tx
			sys.exit()

		sp = SPhoton(self.t0,self.x,self.y,self.z,self.tx,self.ty,self.tz,tau)
		return sp

class Sensor:
	"""
	Describes the sensor (a SiPM)
	"""
	def __init__(self,PDE=0.2,J = 0.1*ns,taur = 0.2*ns,taud = 2*ns):
		self.taud = taud
		self.taur = taur

		self.pde = PDE
		self.J = J

	def PDE(self):
		"""
		The photon detection efficiency
		"""
		return self.pde

	def Jitter(self):
		"""
		Time jitter
		"""
		tj = rnd.gauss(0, self.J*2.355)
		return tj

	def SPE(self,t):
		"""
		Single photoelectron response: t0 is the time of the pe
		"""
		f = exp(-t/self.taud) -exp(-t/self.taur)
		norm = self.taud - self.taur
		return f

	def SensorResponse(self,uvp,STIME):
		"""
		Accumulated response to pe
		"""

		# add the SiPM jitter to the time of the uvp
		t0 = uvp.Time() + self.Jitter()
		
		if DEBUG > 1:
			print "t0 (ps)",t0/ps

		if DEBUG > 2:
			wait()

		if t0 > TAUMAX:  #loop out outliers
			return np.zeros(1)

		sper=[]
		nbins=0

		for tt in STIME:
			t= tt*ps #tt is expressed in ps
			if t > t0:
				spe = self.SPE(t)
			else:
				spe = 0.
			sper.append(spe)

			if DEBUG > 2:
				print "t0, t, spe",t0/ps,tt,spe
			
			nbins+=1

		asper = np.array(sper)
		intg = sum(asper)
		asper = self.PDE()*asper/intg

		if DEBUG > 1:
			print "sum sper ",intg
			print "length of sper vector",len(asper)

		if DEBUG > 2:
			print "sper vector",asper
			wait()

		return asper



class ScintillationEvent:
	"""
	Generates a scintillation event 
	"""
	def __init__(self,WBox,PhotonGenerator,sipm,tau,nphotons,tLevel):
		self.wbox = WBox
		self.pg = PhotonGenerator
		self.np = nphotons
		self.tau = tau
		self.tLevel = tLevel
		self.sipm = sipm
		stime =[]

		for tt in drange(0.,TAUMAX,TBINS):
			stime.append(tt/ps)

		self.STIME = np.array(stime)
		

	def Event(self):
		"""
		The scintillation event 
		"""
		#self.hsper = TH1F("hsper", "spe(au)", NBINS, 0.,(TAUMAX/ps))

		if DEBUG > 0:
			self.hSPER = TH1F("hSPER", "spe(au)", NBINS, 0.,(TAUMAX/ps))
			self.c2 = TCanvas( 'c2', 'SPER', 200, 10, 600, 800 )

		self.SPER = np.zeros(len(self.STIME))

		for i in range(self.np):
		
			if DEBUG > 0:
				print "scintillation photon--> ", i

			uvp = self.pg.GeneratePhoton(tau=self.tau,	
																		x=self.wbox.X()/2,y=self.wbox.Y()/2,
																		z=self.wbox.Z()/2)

			if DEBUG > 0:
				print uvp

			ht0.Fill(uvp.T0()/ns)

			n=0
			while n < MAX_BOUNCES:
				n+=1
				fi = pt.step(uvp)
			
				if DEBUG > 1:
					print "face	index =",fi

				if fi == 0: #hitting x=0 or x=l face, photon is absorbed
					if DEBUG > 1:
						print "hits face 0 or face 2, photon absorbed"

					htime.Fill(uvp.Time()/ns)
					hpath.Fill(uvp.Path()/mm)

					#Efficiency of detection

					test = rnd.uniform(0.,1.)

					if test > DETECTION_EFF:
						if DEBUG> 1:
							print "photon not detected by sensors "

						break

					sper = self.sipm.SensorResponse(uvp,self.STIME)
					if len(sper) == 1: #time too long
						break

					# ibins = min(len(sper),NBINS)
					# for i in xrange(ibins):
					# 	self.hsper.AddBinContent(i, sper[i])

					self.SPER += sper
					break

				else: #is photon absorbed by Teflon?
					
					test = rnd.uniform(0.,1.)
					if test > LXE_REFLECTIVITY:
						if DEBUG > 1:
							print "photon absorbed by Teflon"

						break

					pt.lambert(uvp,fi) 

					if DEBUG > 1:
						print "Lambertian reflection"
			
			hb.Fill(n)
		
		if DEBUG > 0:
			ibins = min(len(self.SPER),NBINS)
			for i in xrange(ibins):
				self.hSPER.AddBinContent(i, self.SPER[i])

		# self.c2.Divide(1,2)
		# self.c2.cd(1)
		# gPad.SetLogy()
		# self.hsper.Draw()
		# self.c2.cd(2)
		# gPad.SetLogy()

		if DEBUG > 0:
			self.hSPER.Draw()
			self.c2.Show()
			wait()
	
	def Trigger(self):
		"""
		Trigger the event 
		"""

		# #print "trigger"
		# ttr = -999
		
		# for ibin in range(NBINS):
		# 	pes = self.hsper.GetBinContent(ibin)
		# 	time = self.hsper.GetBinCenter(ibin)
		# 	#print "time (ps) pes =", time,pes
		# 	if pes > self.tLevel:
		# 		ttr = time
		# 		break

		# print "trigger time1 =",ttr
		# #wait()
		# hCTR1.Fill(ttr)

		ttr = -999
		
		if DEBUG > 1:
			print self.STIME
			print self.SPER

		for ibin in xrange(len(self.STIME)):
			time = self.STIME[ibin]
			pes = self.SPER[ibin]
			if DEBUG > 1:
				print "time (ps) pes =", time,pes

			if pes > self.tLevel:
				ttr = time
				break

		if DEBUG > 0:
			print "trigger time2 ps =",ttr
			wait()
	
		hCTR2.Fill(ttr)


if __name__ == '__main__':

	c1 = TCanvas( 'c1', 'CTR', 200, 10, 600, 800 )
	

	wbox = Box(LX,LY,LZ)
	pt = PhotonTransport(wbox)
	pg = PhotonGenerator(wbox)
	sipm =Sensor(PDE=SIPM_PDE, J=SIPM_J, taur=SIPM_PR,taud=SIPM_PD)
	sce = ScintillationEvent(wbox,pg,sipm,TAU1,NPHOTONS,TLEVEL)

	ht0 = TH1F("ht0", "photon creation time (ns)", NBINS, 0., 3*TAUMAX/ns)
	htime = TH1F("htime", "photon time (ns)", NBINS, 0., 3*TAUMAX/ns)
	hpath = TH1F("hpath", "photon path (mm)", 100, 0., 500.)
	hb = TH1F("hb", "number of bounces", 20, 0., 20.)
	hspe = TH1F("hspe", "spe(au)", NBINS, 0.,3*TAUMAX/ns)
	hCTR1 = TH1F("hCTR1", "ctr1", 200, 800.,1000.)
	hCTR2 = TH1F("hCTR2", "ctr2", 50, 1800.,2200.)
	

	print """
		CTR calculation
		number of events = %d 
		using tau = %7.2f ns
		taumax = %7.2f ns
		number of photons for this tau = %d
		number of time bins for photons = %d
		number of time bins for SPE = %d
		Box dimensions: x = %7.2f mm y = %7.2f mm z = %7.2f mm 
		SiPM PDE = %7.2f
		SiPM Jitter = %7.2f
		SiPM raise constant %7.2f ns
		SiPM decay constant %7.2f ns
		Trigger level = %7.2f
			
	"""%(NEVENTS,TAU/ns,TAUMAX/ns,NG_tau1,NBINS,2*NBINS,
			LX/mm,LY/mm,LZ/mm,SIPM_PDE,SIPM_J,SIPM_PR/ns,SIPM_PD/ns,TLEVEL)

	wait()


	for iev in range(NEVENTS):	
		print "event ", iev
		
		sce.Event()
		sce.Trigger()
		
	c1.Divide(2,2)
	c1.cd(1)
	gPad.SetLogy()
	ht0.Draw()
	
	c1.cd(2)
	gPad.SetLogy()
	htime.Draw()

	c1.cd(3)
	gPad.SetLogy()
	hpath.Draw()

	c1.cd(4)
	gPad.SetLogy(0)
	hCTR2.Draw()

	c1.Show()
	wait()

	# c1.Divide(1,2)
	# c1.cd(1)
	# gPad.SetLogy(0)
	# hCTR1.Draw()
	
	# c1.cd(2)
	# gPad.SetLogy(0)
	# hCTR2.Draw()

	# c1.Show()
	# wait()
	

	#FPULSE.Draw()
	
	# FTAU3.Draw()
	# FTAU2.Draw("same")
	# FTAU1.Draw("same")
	# c1.Show()
	# wait()
	

	
	


