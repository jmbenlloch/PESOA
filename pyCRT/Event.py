from Photon import *
from THisto import *
import random as rnd

log =logging.getLogger("Event")
log.setLevel(logging.INFO)
log.addHandler(ch)
deb = Switch(name="EventDebug",switch=0)
hsw = Switch(name="EventHistos",switch=0)

MAX_BOUNCES=100



class ScintillationEvent:
	"""
	Generates a scintillation event 
	"""
	def __init__(self,lxsc,tauIndex):
		"""
		Scintillation event. it takes:
		the lxsc Box instance 
		Instances to photon generator and transport 
		the interaction point Point3D ixyz  (x0,y0,z0)
		"""
		self.lxsc = lxsc
		self.wbox = lxsc.Box()
		self.lxe = lxsc.LXe()
		self.pg = PhotonGenerator(lxsc)
		self.pt = PhotonTransport(lxsc)
		self.np = (int) (self.lxe.SPhotonsAt511KeV(tauIndex)) #number of photons
		self.sipm = self.lxsc.Sensor()
		self.tauIndex = tauIndex
		
		
		log.info(' Scintillation event: tau = %7.2f ns LXSC box:',
						self.lxe.Lifetime(self.tauIndex)/ns)

		if hsw.Switch() == 1:
			self.th = THisto()
			
			self.th.BookH1("ht0", "photon t0 (ns)", 200, 0., 
							5*self.lxe.Lifetime(tauIndex)/ns)

			self.th.BookH1("htx", "thetax ", 50, -1., 1.)

			self.th.BookH1("hz0", "z0 (mm) ", 50, 0., self.wbox.Z()/mm)
			self.th.BookH1("hz", "z (mm) ", 50, 0., self.wbox.Z()/mm)
			self.th.BookH2("hxy", "xy (mm) ", 50, 0., self.wbox.X()/mm,50, 0., self.wbox.Y()/mm) 
							
			self.th.BookH1("htime", "photon time (ns)", 200, 0., 
							5*self.lxe.Lifetime(tauIndex)/ns)

			self.th.BookH1("hpath", "photon path (mm)", 100, 
							0., 10*self.wbox.X())
			self.th.BookH1("hb", "number of bounces", 50, 0., 50.)
	
		

	def Event(self,ixyz,txyz):
		"""
		The scintillation event
		the interaction point Point3D ixyz  (x0,y0,z0)
		tyxz (tx,ty and tz) are the director cosines of the emited photons 
		"""
		
		PES=[]
		log.debug('On Event: number of photons to shoot = %d ixyz = %s, txyz=%s',
			self.np, ixyz,txyz)

		for i in range(self.np):  # loop over number of photons

			log.debug('scintillation photon--> = %d',i)
			deb.Wait()
			
			#Generate VUV photon
			uvp = self.pg.GeneratePhotonUV(self.tauIndex,
										   x=ixyz.X(),y=ixyz.Y(),z=ixyz.Z(),
										   tx=txyz.X(),ty=txyz.Y(),tz=txyz.Z())


			pe =Photoelectron(uvp)

			if hsw.Switch() == 1:
				self.th.FillH1("ht0",uvp.T0()/ns)
				self.th.FillH1("htx",uvp.DirectionCosines().X())
				self.th.FillH1("hz0",uvp.xyzt0.z/mm)

			log.debug('Generating UV photon = %s',uvp)
			
		
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
							log.debug('photon associated to this pe = %s',
							uvp)
							log.debug('photon associated to this pe: xyzt = %s',
							uvp.xyzt)
							log.debug('photon associated to this pe: xyzt0 = %s',
							uvp.xyzt0)

							deb.Wait()
					
							#histogram path and time
							if hsw.Switch() == 1:
								self.th.FillH1("htime",uvp.Time()/ns)
								self.th.FillH1("hpath",uvp.Path()/mm)
								self.th.FillH2("hxy",uvp.xyzt.x/mm,uvp.xyzt.y/mm)
								self.th.FillH1("hz",uvp.xyzt.z/mm)
							#Sensor Response goes here.
						 
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
							self.th.FillH1("hb",n)
				n+=1
			

		log.info('number of pes = %d ',len(PES))
						
		
		if hsw.Switch() == 1:
			
			self.th.DrawList(["ht0","htime","hz0","hz"],xd=2,yd=2, yscale='log')
			self.th.DrawList(["htx","hxy","hpath","hb"],xd=2,yd=2)
		
		return  PES
	
	

	

	
