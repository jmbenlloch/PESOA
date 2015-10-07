from Photon import *
from ROOT import TH1F, TH2F,TH3F
from ROOT import TCanvas, TPad, TPaveLabel, TPaveText, TColor, TFile
from ROOT import gROOT, gStyle, gPad
gStyle.SetOptStat(1);
gStyle.SetPalette(1);
gStyle.SetCanvasColor(33);
gStyle.SetFrameFillColor(18);
import numpy as np
import random as rnd

log =logging.getLogger("Event")
log.setLevel(logging.INFO)
log.addHandler(ch)
deb = Switch(name="EventDebug",switch=0)
hsw = Switch(name="EventHistos",switch=0)

MAX_BOUNCES=100
		
TBINS = 0.01*ns
TLEVEL = 1 # in pes


class Photoelectron:
	"""
	Represents a photoelectron
	"""
	def __init__(self,phton):
		self.photon = photon

	def Set(self,x,y,z,time):
		"""
		Time and coordinates of the pes.
		For a mask [0,0,1], z= 0 or z = lz, and (x,y) correspond to the sipm coordinates
		"""
		self.x = x
		self.y = y
		self.z = z
		self.time = time
		self.sipm = -1

class ScintillationEvent:
	"""
	Generates a scintillation event 
	"""
	def __init__(self,lxsc,tauIndex,x0,y0,z0):
		"""
		Scintillation event. it takes:
		the lxsc Box instance 
		Instances to photon generator and transport 
		the interaction postion x0,y0,z0 
		"""
		self.lxsc = lxsc
		self.wbox = lxsc.Box()
		self.lxe = lxsc.LXe()
		self.pg = PhotonGenerator(lxsc)
		self.pt = PhotonTransport(lxsc)
		self.np = (int) (self.lxe.SPhotonsAt511KeV(tauIndex))
		self.sipm = self.lxsc.Sensor()
		self.tauIndex = tauIndex
		self.x0 = x0
		self.y0 = y0
		self.z0 = z0
		
		

		log.info(' Scintillation event: tau = %7.2f ns LXSC box: x = y = z = %7.2f mm',
						self.lxe.Lifetime(self.tauIndex)/ns,(self.wbox.X()/2)/mm)

		if hsw.Switch() == 1:
			self.c1 = TCanvas( 'c1', 'CTR', 200, 10, 600, 800 )
			
			self.ht0 = TH1F("ht0", "photon t0 (ns)", 200, 0., 
							5*self.lxe.Lifetime(tauIndex)/ns)

			self.htx = TH1F("htx", "thetax ", 50, -1., 1.) 
							
			self.htime = TH1F("htime", "photon time (ns)", 200, 0., 
							5*self.lxe.Lifetime(tauIndex)/ns)
			self.hpath = TH1F("hpath", "photon path (mm)", 100, 
							0., 10*self.wbox.X())
			self.hb = TH1F("hb", "number of bounces", 50, 0., 50.)
	
		# stime =[]

		# for tt in drange(0.,TAUMAX,TBINS):
		# 	stime.append(tt/ps)

		# self.STIME = np.array(stime)
		

	def Event(self,txp,typ,tzp):
		"""
		The scintillation event
		tx,ty and tz are the director cosines of the emited photons 
		"""
		#self.hsper = TH1F("hsper", "spe(au)", NBINS, 0.,(TAUMAX/ps))

		# if DEBUG > 0:
		# 	self.hSPER = TH1F("hSPER", "spe(au)", NBINS, 0.,(TAUMAX/ps))
		# 	self.c2 = TCanvas( 'c2', 'SPER', 200, 10, 600, 800 )

		# self.SPER = np.zeros(len(self.STIME))


		PES=[]

		for i in range(self.np):  # loop over number of photons

			log.debug('scintillation photon--> = %d',i)
			
			#Generate VUV photon
			uvp = self.pg.GeneratePhotonUV(self.tauIndex,
										   x=self.x0,y=self.y0,z=self.z0,
										   tx=txp,ty=typ,tz=tzp)
			pe =Photoelectron(uvp)

			if hsw.Switch() == 1:
				self.ht0.Fill(uvp.T0()/ns)
				self.htx.Fill(uvp.TX())

			log.debug('Generating UV photon = %s',uvp)
			deb.Wait()

		
			n=0
			while n < MAX_BOUNCES: #propagate up to the max number of bounces
				fi = self.pt.Step(uvp)  #step the photon to the next face
				log.debug('face	index  = %d',fi)

				if self.pt.TestInstrumentedFace(fi) == 1: #hitting instrumented face, photon is absorbed
					log.debug('photon hits instrumented face index  = %d',fi)

					#Efficiency of detection by sensor plane

					if self.pt.TestSensorEfficiency() == True:
						log.debug('photon is detected by sensors')

					#PDE

						if self.pt.TestSensorPDE() == True:
							log.debug('A photoelectron is produced: time stamp = %7.2f ps',
							uvp.Time()/ps)

							deb.Wait()
					
							#histogram path and time
							if hsw.Switch() == 1:
								self.htime.Fill(uvp.Time()/ns)
								self.hpath.Fill(uvp.Path()/mm)
							#Sensor Response goes here.
						
							pe.Set(uvp.X(),uvp.Y(),uvp.Z(),uvp.Time())
							PES.append(pe)

							break

						else:
							log.debug('photon not detected (PDE)')
							break
					else:
						log.debug('photon not detected (Geometry)')
						break

				else: #is VUV photon absorbed by Teflon?
					test = self.pt.TestReflectivityUV()
					
					if test == False:
						log.debug('photon absorbed')
						break
					else: # Lambertian reflection
						log.debug('Lambertian reflection in face index = %d',fi)
						self.pt.Lambert(uvp,fi) 
						if hsw.Switch() == 1:
							self.hb.Fill(n)
				n+=1
			

		log.info('number of pes = %d ',len(PES))
		log.info('time stamp of first photoelectron = %7.2f ps',timeFirstPe/ps)
						
		
		if hsw.Switch() == 1:
			gROOT.FindObject("c1").Update() 
			c1.Divide(2,2)
			c1.cd(1)
			gPad.SetLogy()
			self.ht0.Draw()
	
			c1.cd(2)
			gPad.SetLogy()
			self.htime.Draw()

			c1.cd(3)
			gPad.SetLogy()
			self.hctime.Draw()

			c1.cd(4)
			gPad.SetLogy()
			self.hpath.Draw()

			c1.Show()
			hsw.Wait()

		
		return  PES
	
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

		# ttr = -999
		
		# if DEBUG > 1:
		# 	print self.STIME
		# 	print self.SPER

		# for ibin in xrange(len(self.STIME)):
		# 	time = self.STIME[ibin]
		# 	pes = self.SPER[ibin]
		# 	if DEBUG > 1:
		# 		print "time (ps) pes =", time,pes

		# 	if pes > self.tLevel:
		# 		ttr = time
		# 		break

		# if DEBUG > 0:
		# 	print "trigger time2 ps =",ttr
		# 	wait()
	
		# hCTR2.Fill(ttr)

def GenerateEvent(tauIndex,PDE,NEVENTS):
	"""
	Generates an event 
	"""
	histoPath = "CTRHistos.root"
	fHistoFile = TFile(histoPath,"RECREATE")

	lxe = LXe()
	sipm = Sensor(timeSeries =0, tauMax = 0,
		          pde=PDE,pulseRaise = 1.0*ns, pulseDecay= 3.2*ns,
				  timeJitterFWHM=0*ns)

	plxsc = PLXSC(sigmax=1*mm, sigmay=1*mm, sigmaz=1*mm, sigma0=0.05,
                  sensorEff=0.9, wlsEff=0.8, uvRef=0.95, wlsRef=0.98)

	wbox = Box(50*mm,50*mm,50*mm)
	tpb = WLS()

	imask = [0,0,1]
	lxsc = LXSC(lxe,wbox,tpb,plxsc,sipm,6.2*mm,imask)

	pitch = lxsc.SensorPitch()
	dx = wbox.X()/pitch
	dy = wbox.Y()/pitch 
	nx = (int) (dx)
	ny = (int) (dy) 

	self.xr = self.x + rnd.gauss(0, self.sigmax) # reconstructed x
		self.yr = self.y + rnd.gauss(0, self.sigmay) # reconstructed y
		self.zr = self.z + rnd.gauss(0, self.sigmaz) # reconstructed z

	log.info("pitch = %7.2f mm x = %7.2f mm y = %7.2f mm nx = %d ny = %d",
		pitch,dx,dy,nx,ny)


	hint = TH1F("hint", "frac", 50, 0., wbox.Z()/mm)

	htfe = TH1F("htfe", "time first pe (ps)", 200, 0., 200)
	hxpes = TH1F("hxpes", "x pes in mm", nx, 0., wbox.X()/mm)
	hypes = TH1F("hypes", "y pes in mm", ny, 0., wbox.Y()/mm)
	hxypes = TH2F("hxypes", "x-y pes in mm", nx, 0., wbox.X()/mm,ny, 0., wbox.Y()/mm)
	htpes = TH1F("htpes", "t pes in ns", 200, 0., 5*lxe.Lifetime(tauIndex)/ns)
	
	c1 = TCanvas( 'c1', 'z', 200, 10, 600, 800 )
	

	sce = ScintillationEvent(lxsc,tauIndex)
	

	print """
		CTR calculation
		number of events = %d 
		using tau = %7.2f ns
		number of photons for this tau = %d
		lxsc = %s
		
			
	"""%(NEVENTS,lxe.Lifetime(tauIndex)/ns,sce.np,lxsc)

	wait()


	for iev in range(NEVENTS):	
		z = rnd.uniform(0.,wbox.Z())
		w = lxe.Efficiency(511*keV,z)
		
		log.info('event  = %d: z = %7.2f, w = %7.2f z*w %7.2f',iev, z/mm,w/mm,z*w/mm)
		hint.Fill(z*w/mm)
		#wait()

		PES = sce.Event(wbox.X()/2,wbox.Y()/2,z*w,-9999,-9999,-9999)

		ipes = 0
		timeFirstPe = 1*second

		w = 1./len(PES)
		for pes in PES:
			ipes+=1
			log.debug('pes number  = %d x = %7.2f mm y = %7.2f mm z = %7.2f mm t = %7.2f ps ',
				iev,pes.x,pes.y,pes.z,pes.time)
			photon = pes.photon
			zr = photon.DOI()
			t1 = zr*photon.n/c_light
			t2 = zr*photon.n/c_light

			ctime = photon.time - photon.DOI()*photon.n/c_light
			hxpes.Fill(pes.x/mm,w)
			hypes.Fill(pes.y/mm,w)
			htpes.Fill(pes.time/ns)
			hxypes.Fill(pes.x/mm,pes.y/mm)

			if pes.time < timeFirstPe:
				timeFirstPe = pes.time

		htfe.Fill(timeFirstPe/ps)
			

		log.info('Number of PES, first pes = %d, %7.2f',len(PES), tpe/ps)
		

		#sce.Trigger()
		
	
	c2 = TCanvas( 'c2', 'CTR', 200, 10, 600, 800 )
	c3 = TCanvas( 'c3', 'CTR2', 200, 10, 600, 800 )
	
	c2.Update()
	c2.Divide(2,2)
	c2.cd(1)
	gPad.SetLogy(0)
	hxpes.Draw()
	c2.cd(2)
	gPad.SetLogy(0)
	hypes.Draw()
	c2.cd(3)
	gPad.SetLogy(0)
	hxypes.Draw("Box")
	c2.cd(4)
	gPad.SetLogy(1)
	hint.Draw()
	
	c2.Update()
	wait()

	
	c3.Divide(1,2)
	c3.cd(1)
	gPad.SetLogy(0)
	htfe.Draw()
	c3.cd(2)
	gPad.SetLogy(0)
	htpes.Draw()
	c3.Update()
	
	wait()
	fHistoFile.Write()
  	fHistoFile.Close()
if __name__ == '__main__':
	NEVENTS = 100
	
	GenerateEvent(2,0.2,NEVENTS)
	
	

	
