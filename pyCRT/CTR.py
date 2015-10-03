from math import *
import random as rnd
import numpy as np
#import pylab as P
from LXSC import *
from ROOT import *



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

LXE_REFLECTIVITY = 0.95
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

			uvp = self.pg.GeneratePhoton(tau=self.tau,x=self.wbox.X()/2,y=self.wbox.Y()/2,
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
	

	
	


