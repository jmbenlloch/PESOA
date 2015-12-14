from Centella.AAlgo import AAlgo

from Particles import *
from Util import *
from Geometry import *
from TOF import *
from LXe import *
import random as rnd
import array
from Centella.treeManager import *

"""
This algorithm computes the Coincidence Resolution Time. 
"""


class CRT(AAlgo):

	############################################################
	def __init__(self, param=False, level=1, label="", **kargs):

		"""
		CRT Algorithm
		"""
		#self.m.log(1, 'Constructor()')

		### GENERAL STUFF
		self.name = 'CRT'
		
		AAlgo.__init__(self, param, level, self.name, 0, label, kargs)

    ### PARAMETERS
		
		#self.FTREES = self.strings["FTREES"]
		self.debug = self.ints["Debug"]  #used to stop program at key break points

		#load coordinates of box and fiducial box

		boxCoord1 =self.loadCoord("Box1V")
  		boxCoord2 =self.loadCoord("Box2V")
  		fboxCoord1 =self.loadCoord("FBox1V")
  		fboxCoord2 =self.loadCoord("FBox2V")

  		self.QE = self.doubles["QE"]  #quantum efficiency
  		self.DTMAX = self.doubles["DTMAX"]*ps  #max diff wrt first pes
  		self.SPTR = self.doubles["SPTR"]*ps  #single photon time resolution
  		self.ASIC= self.doubles["ASIC"]*ps  #ASIC contribution
  		 
  		self.NPE = self.ints["NPE"]   #number of pe for time average
  		self.NSIPM = self.ints["NSIPM"]   #number of SiPMs per box

  		#time jitter of SiPM + ASIC
  		self.TJ = sqrt(self.SPTR**2+self.ASIC**2)

  		self.box1ID=self.vints["Box1Id"]  #ids of SiPMs in box1
  		self.box2ID=self.vints["Box2Id"] 

		self.box1 = Box(boxCoord1)
		self.box2 = Box(boxCoord2)
		self.fbox1 = Box(fboxCoord1)
		self.fbox2 = Box(fboxCoord2)
		
		self.m.log(1, "Box1 --", self.box1)
		self.m.log(1, "Box2 ---", self.box2)

		self.m.log(1, "Fiducial Box1 --", self.box1)
		self.m.log(1, "Fiducial Box2 ---", self.box2)
		self.m.log(1, "IDs Box1 --", self.box1ID)
		self.m.log(1, "IDs Box2 --", self.box2ID)

		self.m.log(1, "QE = %7.2f Time Jitter =%7.2f ps  --"%(self.QE, self.TJ/ps))

		self.lxe = LXe() #lxe properties
		print self.lxe

		
		if self.debug == 1:
			wait()


	############################################################		
	def initialize(self):

		self.m.log(1, 'Initialize()')
		self.BookTree()
		
		### Defining histos
		# Event energy histogram

		self.Histos()

    ### Counters:
		self.numInputEvents = 0
		self.numOutputEvents = 0

	#Time Map instance
		self.timeMap = TimeMap(numberOfBoxes = 2, dtmax=self.DTMAX)
		
		return


	############################################################
	def execute(self, event=""):

		self.m.log(2, 'Execute()')		
	
		self.numInputEvents += 1   

		#Fiducial cut. There must be one photon in the fiducial volue
		#(defined by parameters) in each box. return the interaction vertex
		#of the photon

		fiducial = self.Fiducial(event)
		self.hman.fill(self.Fiducial_histo_name,fiducial)
		self.fid[0] = fiducial 

		if self.debug == 1:
			wait()

		if fiducial != 2:
			return False

		
		#Compute a TimeMap including the time-ordered sequence of the
		#PES of each SiPM hit that arrive within DTMAX of the first PE

		self.ComputeTimeMap(event)
		nhitsBox1 = self.timeMap.NumberOfSiPmHits(box=1)
		nhitsBox2 = self.timeMap.NumberOfSiPmHits(box=2)

		self.m.log(2, " nhits box1 = %s nhits box2 = %s"%(
			nhitsBox1,nhitsBox2))

		self.hman.fill(self.nhitsBox1_histo_name,nhitsBox1)
		self.hman.fill(self.nhitsBox2_histo_name,nhitsBox2)
		self.nhitsb1[0] = nhitsBox1
		self.nhitsb2[0] = nhitsBox2

		if self.debug == 1:
			wait()

		if nhitsBox1 < self.NSIPM or nhitsBox2 < self.NSIPM:
			return False
		
		self.m.log(3, " Time maps = %s"%(self.timeMap))
		if self.debug == 1:
			wait()	
		
		#Compute DT 

		dtFirstPe, dtAverageSiPm, dtAverageTime = self.ComputeDT()
		
		self.hman.fill(self.DT_histo_name,dtFirstPe/ps)
		self.hman.fill(self.DTSiPmAvg_histo_name,dtAverageSiPm/ps)
		self.hman.fill(self.DTTimeAvg_histo_name,dtAverageTime/ps)

		self.FirstPeDT[0]=dtFirstPe/ps
		self.AverageSiPmDT[0]=dtAverageSiPm/ps
		self.AverageTimeDT[0]=dtAverageTime/ps
		self.m.log(2,
			"-->dtFirstPe =%7.2f  ps,  "%(dtFirstPe/ps))
		self.m.log(2,
			"-->dtAverageSiPm =%7.2f  ps,  "%(dtAverageSiPm/ps))
		self.m.log(2,
			"-->dtAverageTime =%7.2f  ps,  "%(dtAverageTime/ps))
		
		if self.debug == 1:
			wait()	

		self.numOutputEvents += 1
		self.tman.fill('CRT')
		return True

############################################################
	def finalize(self):

		self.m.log(1, 'Finalize()')

		self.m.log(1, 'Input  Events: ', self.numInputEvents)
		self.m.log(1, 'Output Events: ', self.numOutputEvents)
		

		self.logman["USER"].ints[self.alabel("InputEvents")] = self.numInputEvents
		self.logman["USER"].ints[self.alabel("OutputEvents")] = self.numOutputEvents
	
		#self.tman.save(file_name=self.FTREES)

		return

############################################################
	def Fiducial(self,event):
		"""
		This method checks that there is one photon in the fiducial volue
		(defined by parameters) in each box. It returns 
		fid = 0 if no photon is found in the fiduical volume
		fid = 1 if one photon found in one box
		fid = 2 if one photon found in each box
		vertexBox1, vertexBox2 the vertices of the photons
		"""

		primaryParticles = PrimaryParticles(event)
		self.m.log(2, ' number of primary Particles =%d '%(len(primaryParticles)))
		
		vertexBox1=Point4D()
		vertexBox2=Point4D()
		fid = 0
		for pparticle in primaryParticles:
			self.m.log(2, '\n+++primary particle+++\n')
			ei,ef = particleKineticEnergy(pparticle)

			self.m.log(2,'name = %s t =%7.2f ps E = %7.2f keV'%(
			particleName(pparticle), particleTime(pparticle)/ps,ei/keV))

			x0,y0,z0 = particleInitialVtx(pparticle)
			x,y,z = particleFinalVtx(pparticle)
			
			self.m.log(2, ' x0 =%7.2f mm, y0 =%7.2f mm, z0 =%7.2f mm '%(
			x0/mm,y0/mm,z0/mm))
			self.m.log(2, ' xf =%7.2f mm, yf =%7.2f mm, zf =%7.2f mm '%(
			x/mm,y/mm,z/mm))

			if self.fbox1.Active((x,y,z)) == True:
				self.m.log(2,'gamma found in box1')
				
				self.hman.fill(self.XYZBox1_histo_name, 
					x/mm,y/mm,z/mm)
				self.hman.fill(self.XYBox1_histo_name, 
					x/mm,y/mm)
				self.hman.fill(self.ZBox1_histo_name, 
					z/mm)

				vertexBox1.x = x
				vertexBox1.y = y
				vertexBox1.z = z
				vertexBox1.t = particleTime(pparticle)

				self.xb1[0]=x/mm
				self.yb1[0]=y/mm
				self.zb1[0]=z/mm
				self.tpb1[0]=vertexBox1.t/ps

				fid+=1
			
			elif self.fbox2.Active((x,y,z)) == True:
				self.m.log(2,'gamma found in box2')
				
				self.hman.fill(self.XYZBox2_histo_name, 
					x/mm,y/mm,z/mm)
				self.hman.fill(self.XYBox2_histo_name, 
					x/mm,y/mm)
				self.hman.fill(self.ZBox2_histo_name, 
					z/mm)

				vertexBox2.x = x
				vertexBox2.y = y
				vertexBox2.z = z
				vertexBox2.t = particleTime(pparticle)

				self.xb2[0]=x/mm
				self.yb2[0]=y/mm
				self.zb2[0]=z/mm
				self.tpb2[0]=vertexBox2.t/ps

				fid+=1
			else:
				self.m.log(2,'gamma not found in box1 or in box2')
				self.m.log(2, ' x0 =%7.2f mm, y0 =%7.2f mm, z0 =%7.2f mm '%(
			x0/mm,y0/mm,z0/mm))
				self.m.log(2, ' xf =%7.2f mm, yf =%7.2f mm, zf =%7.2f mm '%(
			x/mm,y/mm,z/mm))
				self.m.log(2, "fBox1 --", self.fbox1)
				self.m.log(2, "fBox2 ---", self.fbox2)
				
		self.timeMap.SetInteractionVertices((vertexBox1,vertexBox2))
		return fid

############################################################
	def ComputeTimeMap(self,event):
		"""
		Compute time maps
		"""

		sensorhits =  event.GetMCSensHits()
		self.m.log(3, " Compute time map: event has %d sensor hits = "%(
			len(sensorhits)))
		

		timeMapBox1={}
		timeMapBox2={}

		ng1 = 0
		ng2 = 0
		for hit in sensorhits:
			hid = hit.GetSensorID()
			xh = hit.GetPosition().x()
			yh = hit.GetPosition().y()
			zh = hit.GetPosition().z()
			Ah = hit.GetAmplitude()

			# time first pes (MC pes, QE =1)
			waveform = hit.GetWaveform().GetData()

			self.m.log(5, " waveform for hit ID = %d: length = %d"%(
				hid,len(waveform)))
			
			timeBin0 = waveform[0]
			tbin = timeBin0.first
			time0 = tbin*5*ps

			sipmhit = SiPMHit(hid,xh,yh,zh,Ah,time0,self.QE,self.DTMAX)

			self.m.log(5, " ++++++hit = %s"%(sipmhit))
			
			if self.debug == 2:
				wait()
			
			if self.boxId(hid) == 1:
				
				ng1+=Ah
				self.photSiPMb1[0]=Ah
			else:
				ng2+=Ah
				self.photSiPMb2[0]=Ah
			
			# keep all the pes within DTMAX (~200 ps) of first pe

			np=0
			Q=0
			for timeBins in waveform:
				# get the arrival time of the pe
				tbin = timeBins.first
				time = tbin*5*ps
				A = timeBins.second
				np+=1
				
				self.m.log(6, "pe number = %d tbin = %d time = %7.2f ps A = %7.2f pes"%(
					np, tbin,time/ps,A))

				self.m.log(6, " DT wrt 1st = %7.2f ps"%((time - time0)/ps))
					
				#A is the number of pes in the bin. It is equal to 1
				# most of the time (because the time bining id very thin)
				#but it can be a higher number. A is computed assuming QE=1
				# One needs to correct for the QE of the device. 

				npes = 0
				for p in xrange(0,A):
					if rnd.uniform(0.,1.) <= self.QE:
						npes+=1
						Q+=npes

				self.m.log(6, " npes = %d, Q = %d "%(
					npes, Q))

				#is there a real pe in this time bin? otherwise loop
				if npes < 1:
					continue  #didn't pass because of finite QE
				else: #valid pes, keep the time stamp if time within DTMAX

					if abs(time - time0) < self.DTMAX or len(sipmhit.W)==0:
						sipmhit.W.append(time)
					else:
						break 
			
			#sipmhit.Q = Q # one can run over the waveform and get Q
			sipmhit.Q = Ah*self.QE	
			self.m.log(5, " ++++++hit again = %s"%(sipmhit))

			if Q == 0: # QE killed this boy
				continue

			if self.debug == 2:
				wait()

			if self.boxId(hid) == 1:
				timeMapBox1[hid] = sipmhit
				self.pesSiPMb1[0] = sipmhit.NumberOfPE()
				self.dtPESb1[0] = (time - time0)/ps
				self.hman.fill(self.npesBOX1_histo_name,sipmhit.NumberOfPE())

			elif self.boxId(hid) == 2:
				timeMapBox2[hid] = sipmhit 
				self.pesSiPMb2[0] = sipmhit.NumberOfPE()
				self.dtPESb2[0] = (time - time0)/ps
				self.hman.fill(self.npesBOX2_histo_name,sipmhit.NumberOfPE())
			else:
				print "error: hit index not in range, hid=%d"%(hid)
				sys.exit(0)
			

		self.hman.fill(self.NGBOX1_histo_name,ng1)
		self.hman.fill(self.NGBOX2_histo_name,ng2)
		self.nPhotb1[0]=ng1
		self.nPhotb2[0]=ng2 
		
		if ng1 ==0 or ng2 ==0:
			print "!!! ng1 = %d, ng2=%d"%(ng1,ng2)
			return False

		# sort the maps according to the time stamp of first pe
		TimeMapBox1 = sorted(timeMapBox1.items(), key=sortSiPmHits)
		TimeMapBox2 = sorted(timeMapBox2.items(), key=sortSiPmHits)

		self.m.log(3, ' ng1 =%7.2f, ng2 =%7.2f '%(ng1,ng2))
		self.m.log(5, " time map box1 = ",timeMapBox1)
		self.m.log(5, " time map box2 = ",timeMapBox2)

		self.timeMap.SetSiPmMaps((TimeMapBox1,TimeMapBox2))
		
		self.m.log(5, " event TimeMap = %s"%(self.timeMap))
		
		if self.debug == 1:
			wait()

###########################################################
	def ComputeDT(self):
		"""
		Compute DT using three methods:
		1) FirstPE to arrive to SiPMs
		2) Average of the first pes to arrive to any siPM
		3) Average of the first wave of pes to arrive to any Sipm 
		"""

		self.m.log(3," ---ComputeDT---")

		vertexBox1 =self.timeMap.InteractionVertex(box=1)
		vertexBox2 =self.timeMap.InteractionVertex(box=2)

		DT1 = self.ComputeDT1((0,0,0),vertexBox1.XYZ(),vertexBox2.XYZ())
		self.m.log(3," DT1 =  %7.2f ps "%(DT1/ps))

		self.DT1[0]=DT1/ps
		self.hman.fill(self.DT1_histo_name,DT1/ps)

		DT2,dtHit12 = self.DTFirstPe()

		dtFirstPe  = self.DTOF(dtHit12,DT1,DT2)	
		self.m.log(3,
			"dtFirstPe =%7.2f  ps,  "%(dtFirstPe/ps))

		DT2,dtHit12 = self.DTAverageSiPm()
		dtAverageSiPm  = self.DTOF(dtHit12,DT1,DT2)	
		self.m.log(3,
			"dtAverageSiPm =%7.2f  ps,  "%(dtAverageSiPm/ps))
		
		DT2,dtHit12 = self.DTAverageTime()
		dtAverageTime = self.DTOF(dtHit12,DT1,DT2)
		self.m.log(3,
			"dtAverageTime =%7.2f  ps,  "%(dtAverageTime/ps))

		return dtFirstPe,dtAverageSiPm,dtAverageTime

		
###########################################################
	def DTOF(self,dthit12,dt1,dt2):
		"""
		Compute difference of time of flight: DTOF
		"""
		return (dthit12 - dt1 - dt2)/2.

###########################################################
	def DTFirstPe(self):
		"""
		Compute DT using the first PE
		"""

		self.m.log(3," ---DTFirstPe---")

		siPMHit1 = self.timeMap.SiPmHit(box=1,index=0)
		siPMHit2 = self.timeMap.SiPmHit(box=2,index=0)
		vertexBox1 =self.timeMap.InteractionVertex(box=1)
		vertexBox2 =self.timeMap.InteractionVertex(box=2)

		self.m.log(3," DTFirstPe: Hit 1 =  %s"%(siPMHit1))
		self.m.log(3," DTFirstPe: Hit 2 =  %s"%(siPMHit2))

		DT2 = self.ComputeDT2(siPMHit1.XYZ(),vertexBox1.XYZ(),
								   siPMHit2.XYZ(),vertexBox2.XYZ())
		self.m.log(3," DT2 =  %7.2f ps "%(DT2/ps))

		dtHit12 = self.ComputeDTHit12(siPMHit1.TimeFirstPE(), siPMHit2.TimeFirstPE())  
		self.m.log(3,"dtHit12 =%7.2f  ps,  "%(dtHit12/ps))

			
		self.tFstSiPMb1[0] = siPMHit1.TimeFirstPE()/ps
		self.tFstSiPMb2[0] = siPMHit2.TimeFirstPE()/ps			
		self.DT2[0]=DT2/ps
		self.DTHit12[0]=dtHit12/ps

		self.hman.fill(self.T0Box1_histo_name,siPMHit1.TimeFirstPE()/ps)		
		self.hman.fill(self.DT2_histo_name,DT2/ps)
		self.hman.fill(self.DTHit12_histo_name,dtHit12/ps)

		return DT2,dtHit12

###########################################################
	def DTAverageSiPm(self):
		"""
		Compute DT using the first PE of all the SiPMs with a first PE
		within DTMAX of the first PE. 
		"""

		self.m.log(3," ---DTAverageSiPm---") 

		sipmDTMap1 = self.timeMap.SiPmMapDT(box=1)
		sipmDTMap2 = self.timeMap.SiPmMapDT(box=2)
		vertexBox1 =self.timeMap.InteractionVertex(box=1)
		vertexBox2 =self.timeMap.InteractionVertex(box=2)
	
		ll = min(len(sipmDTMap1),len(sipmDTMap2))

		dt2 = 0
		dth12 = 0
		for idx in xrange(0,ll):
			siPMHit1 = self.timeMap.SiPmHit(box=1,index=idx)
			siPMHit2 = self.timeMap.SiPmHit(box=2,index=idx)
			dt2+= self.ComputeDT2(siPMHit1.XYZ(),vertexBox1.XYZ(),
								   siPMHit2.XYZ(),vertexBox2.XYZ())

			dth12+=self.ComputeDTHit12(siPMHit1.TimeFirstPE(), siPMHit2.TimeFirstPE())

		DT2 = dt2/float(ll)
		dtHit12 =dth12/float(ll)

		self.m.log(3," DTSiPmAvg2 =  %7.2f ps "%(DT2/ps)) 
		self.m.log(3,"dtSiPmAvgHit12 =%7.2f  ps,  "%(dtHit12/ps))

		self.DTSiPmAvg2[0]=DT2/ps
		self.DTSiPmAvgHit12[0]=dtHit12/ps
		
		self.hman.fill(self.DTSiPmAvg2_histo_name,DT2/ps)
		self.hman.fill(self.DTSiPmAvgHit12_histo_name,dtHit12/ps)


		return DT2,dtHit12


###########################################################
	def DTAverageTime(self):
		"""
		Compute DT using all the PEs in the event with time stamp 
		less than DTMAX. 
		"""

		self.m.log(3," ---DTAverageTime---") 
		
		sipmTimeMap1 = self.timeMap.SiPmTimeMap(box=1)
		sipmTimeMap2 = self.timeMap.SiPmTimeMap(box=2)
		vertexBox1 =self.timeMap.InteractionVertex(box=1)
		vertexBox2 =self.timeMap.InteractionVertex(box=2)
		
		ll = min(len(sipmTimeMap1),len(sipmTimeMap2))

		dt2 = 0
		dth12 = 0
		for idx in xrange(0,ll):
			sipmTime1 = sipmTimeMap1[idx][0]
			siPMHit1 = sipmTimeMap1[idx][1]
			sipmTime2 = sipmTimeMap2[idx][0]
			siPMHit2 = sipmTimeMap2[idx][1]
			
			dt2+= self.ComputeDT2(siPMHit1.XYZ(),vertexBox1.XYZ(),
								   siPMHit2.XYZ(),vertexBox2.XYZ())

			dth12+=self.ComputeDTHit12(sipmTime1,sipmTime2)

		DT2 = dt2/float(ll)
		dtHit12 =dth12/float(ll)

		self.m.log(3," DTTimeAvg2 =  %7.2f ps "%(DT2/ps)) 
		self.m.log(3,"dtTimeAvgHit12 =%7.2f  ps,  "%(dtHit12/ps))

		self.DTTimeAvg2[0]=DT2/ps
		self.DTTimeAvgHit12[0]=dtHit12/ps
		
		self.hman.fill(self.DTTimeAvg2_histo_name,DT2/ps)
		self.hman.fill(self.DTTimeAvgHit12_histo_name,dtHit12/ps)

		return DT2,dtHit12

###########################################################
	def ComputeDT1(self,v0,v1,v2):
		"""
		Compute DT1 --> difference of TOF for the two gammas
		between vertex and interaction point 
		"""

		d0v1 = distance(v0,v1)
		d0v2 = distance(v0,v2)
		t0v1 = d0v1/c_light
		t0v2 = d0v2/c_light
		DT1 = t0v1 -t0v2

		self.m.log(4," d0v1 = %7.2f mm, t0v1 = %7.2f ps "%(d0v1/mm,t0v1/ps))
		self.m.log(4," d0v2 = %7.2f mm, t0v2 = %7.2f ps "%(d0v2/mm,t0v2/ps))

		return DT1

###########################################################
	def ComputeDT2(self,vhit1,vertex1,vhit2,vertex2):
		"""
		Compute DT2 --> difference of TOF for the VUV photons 
		propagating from interaction vertex to SiPM in both boxes
		"""


		dbox1 = distance(vhit1,vertex1)
		tpath1 = dbox1*self.lxe.RefractionIndexUV()/c_light
		dbox2 = distance(vhit2,vertex2)
		tpath2 = dbox2*self.lxe.RefractionIndexUV()/c_light
		DT2 = tpath1 - tpath2

		self.hman.fill(self.TBox1_histo_name,tpath1/ps)
		
		self.m.log(4,
			"dbox1 =%7.2f mm,tpath1= %7.2f ps,  "%(dbox1/mm,tpath1/ps))
		self.m.log(4,
			"dbox2 =%7.2f mm,tpath2= %7.2f ps,  "%(dbox2/mm,tpath2/ps))

		return DT2
###########################################################
	def ComputeDTHit12(self,thitBox1, thitBox2):
		"""
		Compute DTHit12: time difference between hit in box1 
		and hit in box2
		"""

		dtHit12 = thitBox1 - thitBox2
		return dtHit12


	############################################################
	def particleInfo(self, lvl, particle):
		self.m.log(lvl,' ++ParticleInfo++')
		
		ti,tf = particleKineticEnergy(particle)
		self.m.log(lvl,'Ti = %7.2f keV Tf= %7.2f keV '%(
			ti/keV,tf/keV)) 

		x,y,z = particleInitialVtx(particle)
		self.m.log(lvl, ' Initial vertex: x =%7.2f mm, y =%7.2f mm, z =%7.2f mm '%(
			x/mm,y/mm,z/mm))
		x,y,z = particleFinalVtx(particle)
		self.m.log(lvl, ' Final vertex: x =%7.2f mm, y =%7.2f mm, z =%7.2f mm '%(
			x/mm,y/mm,z/mm))

		px,py,pz = particleInitialMomentum(particle)
		self.m.log(lvl, ' Initial momentum: px =%7.2f keV, py =%7.2f keV, z =%7.2f keV '%(
			px/keV,py/keV,pz/keV))

		px,py,pz = particleFinalMomentum(particle)
		self.m.log(lvl, ' Final momentum: px =%7.2f keV, py =%7.2f keV, z =%7.2f keV '%(
			px/keV,py/keV,pz/keV))

	

	############################################################
	def loadCoord(self, blb):
		boxCoord =[]
		for i in range(1,9):
			lbl = blb+str(i)
			self.m.log(3, "loading parameter", lbl)
			coord  = self.vdoubles[lbl]
			try:
				boxCoord.append(coord);
			except KeyError:
				self.m.log(1, "WARNING!! Parameter %s not defined."%(lbl))
				exit(0)
		return boxCoord

	############################################################
	def boxId(self, hid):
		b1f1min = self.box1ID[0]
		b1f1max = self.box1ID[1]
		b1f2min = self.box1ID[2]
		b1f2max = self.box1ID[3]

		b2f1min = self.box2ID[0]
		b2f1max = self.box2ID[1]
		b2f2min = self.box2ID[2]
		b2f2max = self.box2ID[3]

		if inRange(hid,b1f1min,b1f1max) or inRange(hid,b1f2min,b1f2max):
			return 1

		if inRange(hid,b2f1min,b2f1max) or inRange(hid,b2f2min,b2f2max):
			return 2

		return 0

############################################################
	def BookTree(self):
	
		self.fid = array.array('f',[0.]) #0,1 or 2 gammas found in boxes

		self.xb1= array.array('f',[0.])   # (x,y,z) vertex in box1
		self.yb1= array.array('f',[0.])
		self.zb1= array.array('f',[0.])
		self.tpb1= array.array('f',[0.])  # time of particle in box1
		self.photSiPMb1= array.array('f',[0.])  # photons/SiPM
		self.pesSiPMb1= array.array('f',[0.])  # pes/SiPM
		self.nPhotb1= array.array('f',[0.])  # total photons in b1
		self.nhitsb1 = array.array('f',[0.]) #number of SiPM with signal in box1
		self.tFstSiPMb1 = array.array('f',[0.]) #time of first SiPM in box1
		self.dtPESb1  = array.array('f',[0.]) #dt of pes with first SiPM in box1
		
		
		self.xb2= array.array('f',[0.])   # (x,y,z) vertex in box2
		self.yb2= array.array('f',[0.])
		self.zb2= array.array('f',[0.])
		self.tpb2= array.array('f',[0.])  # time of particle in box2
		self.photSiPMb2= array.array('f',[0.])
		self.pesSiPMb2= array.array('f',[0.])  # photons/SiPM
		self.nPhotb2= array.array('f',[0.])
		self.nhitsb2 = array.array('f',[0.])
		self.tFstSiPMb2 = array.array('f',[0.])
		self.dtPESb2  = array.array('f',[0.])
		
		self.DT1 = array.array('f',[0.]) # (d(0,v1) - d(0,v2))/c

		self.DT2 = array.array('f',[0.]) #(d(v1,hit1) - d(v2,hit2))*(n/c)
		self.DTHit12 = array.array('f',[0.]) # time (hit1 - hit2)
		self.DTSiPmAvg2 = array.array('f',[0.]) 
		self.DTSiPmAvgHit12 = array.array('f',[0.]) 
		self.DTTimeAvg2 = array.array('f',[0.]) 
		self.DTTimeAvgHit12 = array.array('f',[0.]) 
		
		self.FirstPeDT =array.array('f',[0.])
		self.AverageSiPmDT =array.array('f',[0.])
		self.AverageTimeDT =array.array('f',[0.])
		
		self.tman.book('CRT',"Coincidence Resoltion Time")
		self.tman.addBranch('CRT','fid',self.fid,dim=1)
		self.tman.addBranch('CRT','DT1',self.DT1,dim=1)

		self.tman.addBranch('CRT','DT2',self.DT2,dim=1)
		self.tman.addBranch('CRT','DTHit12',self.DTHit12,dim=1)
		self.tman.addBranch('CRT','DTSiPmAvg2',self.DTSiPmAvg2,dim=1)
		self.tman.addBranch('CRT','DTSiPmAvgHit12',self.DTSiPmAvgHit12,dim=1)
		self.tman.addBranch('CRT','DTTimeAvg2',self.DTTimeAvg2,dim=1)
		self.tman.addBranch('CRT','DTTimeAvgHit12',self.DTTimeAvgHit12,dim=1)

		self.tman.addBranch('CRT','DTFirstPe',self.FirstPeDT,dim=1)
		self.tman.addBranch('CRT','DTAverageSiPm',self.AverageSiPmDT,dim=1)
		self.tman.addBranch('CRT','DTAverageTime',self.AverageTimeDT,dim=1)

		self.tman.addBranch('CRT','xb1',self.xb1,dim=1)
		self.tman.addBranch('CRT','yb1',self.yb1,dim=1)
		self.tman.addBranch('CRT','zb1',self.zb1,dim=1)
		self.tman.addBranch('CRT','tGammab1',self.tpb1,dim=1)
		self.tman.addBranch('CRT','photSiPMb1',self.photSiPMb1,dim=1)
		self.tman.addBranch('CRT','pesSiPMb1',self.pesSiPMb1,dim=1)
		self.tman.addBranch('CRT','tFstSiPMb1',self.tFstSiPMb1,dim=1)
		self.tman.addBranch('CRT','dtPESb1',self.dtPESb1,dim=1)

		self.tman.addBranch('CRT','nPhotb1',self.nPhotb1,dim=1)
		self.tman.addBranch('CRT','nhitsb1',self.nhitsb1,dim=1)

		self.tman.addBranch('CRT','xb2',self.xb2,dim=1)
		self.tman.addBranch('CRT','yb2',self.yb2,dim=1)
		self.tman.addBranch('CRT','zb2',self.zb2,dim=1)
		self.tman.addBranch('CRT','tGammab2',self.tpb2,dim=1)
		self.tman.addBranch('CRT','photSiPMb2',self.photSiPMb2,dim=1)
		self.tman.addBranch('CRT','pesSiPMb2',self.pesSiPMb2,dim=1)
		self.tman.addBranch('CRT','tFstSiPMb2',self.tFstSiPMb2,dim=1)
		self.tman.addBranch('CRT','dtPESb2',self.dtPESb2,dim=1)
		self.tman.addBranch('CRT','nPhotb2',self.nPhotb2,dim=1)
		self.tman.addBranch('CRT','nhitsb2',self.nhitsb2,dim=1)
	
   
	############################################################		
	def Histos(self):
		"""
		book the histograms for the algo 
		"""
		Fiducial_histo_desc = "Fiducial"
		self.Fiducial_histo_name = self.alabel(Fiducial_histo_desc)
		self.hman.h1(self.Fiducial_histo_name, Fiducial_histo_desc, 
			10, 0, 5)
		self.hman.fetch(
			self.Fiducial_histo_name).GetXaxis().SetTitle(
			"Fiducial interactions")

		NGBOX1_histo_desc = "NGBOX1"
		self.NGBOX1_histo_name = self.alabel(NGBOX1_histo_desc)
		self.hman.h1(self.NGBOX1_histo_name, NGBOX1_histo_desc, 
			100, 0, 40000)
		self.hman.fetch(
			self.NGBOX1_histo_name).GetXaxis().SetTitle(
			"Number of gammas in box1")

		NGBOX2_histo_desc = "NGBOX2"
		self.NGBOX2_histo_name = self.alabel(NGBOX2_histo_desc)
		self.hman.h1(self.NGBOX2_histo_name, NGBOX2_histo_desc, 
			100, 0, 40000)
		self.hman.fetch(
			self.NGBOX2_histo_name).GetXaxis().SetTitle(
			"Number of gammas in box2")

		npesBOX1_histo_desc = "NPESBOX1"
		self.npesBOX1_histo_name = self.alabel(npesBOX1_histo_desc)
		self.hman.h1(self.npesBOX1_histo_name, npesBOX1_histo_desc, 
			100, 0, 500)
		self.hman.fetch(
			self.npesBOX1_histo_name).GetXaxis().SetTitle(
			"Number of pes/SiPM in box1")

		npesBOX2_histo_desc = "NPESBOX2"
		self.npesBOX2_histo_name = self.alabel(npesBOX2_histo_desc)
		self.hman.h1(self.npesBOX2_histo_name, npesBOX2_histo_desc, 
			100, 0, 500)
		self.hman.fetch(
			self.npesBOX2_histo_name).GetXaxis().SetTitle(
			"Number of pes/SiPM in box2")

		nhitsBox1_histo_desc = "nhitsBox1"
		self.nhitsBox1_histo_name = self.alabel(nhitsBox1_histo_desc)
		self.hman.h1(self.nhitsBox1_histo_name, nhitsBox1_histo_desc, 
			100, 0, 200)
		self.hman.fetch(
			self.nhitsBox1_histo_name).GetXaxis().SetTitle(
			"Number of sensor hits in box1")

		nhitsBox2_histo_desc = "nhitsBox2"
		self.nhitsBox2_histo_name = self.alabel(nhitsBox2_histo_desc)
		self.hman.h1(self.nhitsBox2_histo_name, nhitsBox2_histo_desc, 
			100, 0, 200)
		self.hman.fetch(
			self.nhitsBox2_histo_name).GetXaxis().SetTitle(
			"Number of sensor hits in box2")


		XYZBox1_histo_desc = "XYZBox1"
		self.XYZBox1_histo_name = self.alabel(XYZBox1_histo_desc)
		self.hman.h3(self.XYZBox1_histo_name, XYZBox1_histo_desc, 
			10, self.box1.xmin, self.box1.xmax,
           	10, self.box1.ymin, self.box1.ymax,
           	10, self.box1.zmin, self.box1.zmax)
		self.hman.fetch(
			self.XYZBox1_histo_name).GetXaxis().SetTitle(
			"XYZ interaction box 1 (mm)")

		XYBox1_histo_desc = "XYBox1"
		self.XYBox1_histo_name = self.alabel(XYBox1_histo_desc)
		self.hman.h2(self.XYBox1_histo_name, XYBox1_histo_desc, 
			25, self.box1.xmin, self.box1.xmax,
           	25, self.box1.ymin, self.box1.ymax)  	
		self.hman.fetch(
			self.XYBox1_histo_name).GetXaxis().SetTitle(
			"XY interaction box 1 (mm)")

		ZBox1_histo_desc = "ZBox1"
		self.ZBox1_histo_name = self.alabel(ZBox1_histo_desc)
		self.hman.h1(self.ZBox1_histo_name, ZBox1_histo_desc, 
			25, self.box1.zmin, self.box1.zmax)
		self.hman.fetch(
			self.ZBox1_histo_name).GetXaxis().SetTitle(
			"Z interaction box 1 (mm)")

		
		XYZBox2_histo_desc = "XYZBox2"
		self.XYZBox2_histo_name = self.alabel(XYZBox2_histo_desc)
		self.hman.h3(self.XYZBox2_histo_name, XYZBox2_histo_desc, 
			10, self.box2.xmin, self.box2.xmax,
           	10, self.box2.ymin, self.box2.ymax,
           	10, self.box2.zmin, self.box2.zmax)
		self.hman.fetch(
			self.XYZBox2_histo_name).GetXaxis().SetTitle(
			"XYZ interaction box 2 (mm)")

		XYBox2_histo_desc = "XYBox2"
		self.XYBox2_histo_name = self.alabel(XYBox2_histo_desc)
		self.hman.h2(self.XYBox2_histo_name, XYBox2_histo_desc, 
			25, self.box2.xmin, self.box2.xmax,
           	25, self.box2.ymin, self.box2.ymax)
           	
		self.hman.fetch(
			self.XYBox2_histo_name).GetXaxis().SetTitle(
			"XY interaction box 2 (mm)")

		ZBox2_histo_desc = "ZBox2"
		self.ZBox2_histo_name = self.alabel(ZBox2_histo_desc)
		self.hman.h1(self.ZBox2_histo_name, ZBox2_histo_desc, 
			25, self.box1.zmin, self.box1.zmax)
   
		self.hman.fetch(
			self.ZBox2_histo_name).GetXaxis().SetTitle(
			"Z interaction box 2 (mm)")

		T0Box1_histo_desc = "T0Box1"
		self.T0Box1_histo_name = self.alabel(T0Box1_histo_desc)
		self.hman.h1(self.T0Box1_histo_name, T0Box1_histo_desc, 
			50, 300, 800)
		self.hman.fetch(
			self.T0Box1_histo_name).GetXaxis().SetTitle(
			"t hit box1 (ps)")

		TBox1_histo_desc = "TBox1"
		self.TBox1_histo_name = self.alabel(TBox1_histo_desc)
		self.hman.h1(self.TBox1_histo_name, TBox1_histo_desc, 
			50, 0, 500)
		self.hman.fetch(
			self.TBox1_histo_name).GetXaxis().SetTitle(
			"time T0 hit to vertex box 1 (ps)")


		DT1_histo_desc = "DT1"
		self.DT1_histo_name = self.alabel(DT1_histo_desc)
		self.hman.h1(self.DT1_histo_name, DT1_histo_desc, 
			50, -500, 500)
		self.hman.fetch(
			self.DT1_histo_name).GetXaxis().SetTitle(
			"(d(0,v1) - d(0,v2))/c (ps)")

		DT2_histo_desc = "DT2"
		self.DT2_histo_name = self.alabel(DT2_histo_desc)
		self.hman.h1(self.DT2_histo_name, DT2_histo_desc, 
			50, -500, 500)
		self.hman.fetch(
			self.DT2_histo_name).GetXaxis().SetTitle(
			"(d(v1-h1) - d(v2 -h2))*(n/c) (ps)")

		DTSiPmAvg2_histo_desc = "DTSiPmAvg2"
		self.DTSiPmAvg2_histo_name = self.alabel(DTSiPmAvg2_histo_desc)
		self.hman.h1(self.DTSiPmAvg2_histo_name, DTSiPmAvg2_histo_desc, 
			50, -500, 500)
		self.hman.fetch(
			self.DTSiPmAvg2_histo_name).GetXaxis().SetTitle(
			"(d(v1-h1) - d(v2 -h2))*(n/c) (ps)")

		DTTimeAvg2_histo_desc = "DTTimeAvg2"
		self.DTTimeAvg2_histo_name = self.alabel(DTTimeAvg2_histo_desc)
		self.hman.h1(self.DTTimeAvg2_histo_name, DTTimeAvg2_histo_desc, 
			50, -500, 500)
		self.hman.fetch(
			self.DTTimeAvg2_histo_name).GetXaxis().SetTitle(
			"(d(v1-h1) - d(v2 -h2))*(n/c) (ps)")

		DTHit12_histo_desc = "DTHit12"
		self.DTHit12_histo_name = self.alabel(DTHit12_histo_desc)
		self.hman.h1(self.DTHit12_histo_name, DTHit12_histo_desc, 
			50, -500, 500)
		self.hman.fetch(
			self.DTHit12_histo_name).GetXaxis().SetTitle(
			"( thit1 - thit2 (ps)")

		DTSiPmAvgHit12_histo_desc = "DTSiPmAvgHit12"
		self.DTSiPmAvgHit12_histo_name = self.alabel(DTSiPmAvgHit12_histo_desc)
		self.hman.h1(self.DTSiPmAvgHit12_histo_name, DTSiPmAvgHit12_histo_desc, 
			50, -500, 500)
		self.hman.fetch(
			self.DTSiPmAvgHit12_histo_name).GetXaxis().SetTitle(
			"( thit1 - thit2 (ps)")

		DTTimeAvgHit12_histo_desc = "DTTimeAvgHit12"
		self.DTTimeAvgHit12_histo_name = self.alabel(DTTimeAvgHit12_histo_desc)
		self.hman.h1(self.DTTimeAvgHit12_histo_name, DTTimeAvgHit12_histo_desc, 
			50, -500, 500)
		self.hman.fetch(
			self.DTTimeAvgHit12_histo_name).GetXaxis().SetTitle(
			"( thit1 - thit2 (ps)")

		DT_histo_desc = "DTFirstPe"
		self.DT_histo_name = self.alabel(DT_histo_desc)
		self.hman.h1(self.DT_histo_name, DT_histo_desc, 
			50, -200, 200)
		self.hman.fetch(
			self.DT_histo_name).GetXaxis().SetTitle(
			"DT First Pe (ps)")

		DTSiPmAvg_histo_desc = "DTSiPmAvg"
		self.DTSiPmAvg_histo_name = self.alabel(DTSiPmAvg_histo_desc)
		self.hman.h1(self.DTSiPmAvg_histo_name, DTSiPmAvg_histo_desc, 
			50, -200, 200)
		self.hman.fetch(
			self.DTSiPmAvg_histo_name).GetXaxis().SetTitle(
			"DTSiPmAvg (ps)")

		DTTimeAvg_histo_desc = "DTTimeAvg"
		self.DTTimeAvg_histo_name = self.alabel(DTTimeAvg_histo_desc)
		self.hman.h1(self.DTTimeAvg_histo_name, DTTimeAvg_histo_desc, 
			50, -200, 200)
		self.hman.fetch(
			self.DTTimeAvg_histo_name).GetXaxis().SetTitle(
			"DTTimeAvg (ps)")

		

	
	