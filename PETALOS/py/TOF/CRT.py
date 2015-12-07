from Centella.AAlgo import AAlgo

from Particles import *
from Util import *
from Geometry import *
from TOF import *
from LXe import *
import random as rnd

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
    # Box coordinates

  		#print self.vdoubles
  		self.debug = self.ints["Debug"]
  		boxCoord1 =self.loadCoord("Box1V")
  		boxCoord2 =self.loadCoord("Box2V")
  		fboxCoord1 =self.loadCoord("FBox1V")
  		fboxCoord2 =self.loadCoord("FBox2V")
  		self.QE = self.doubles["QE"]  #quantum efficiency
  		self.TJ = self.doubles["JITTER"]*ps  #time jitter of SiPM + ASIC
  		self.NPE = self.ints["NPE"]   #number of pe for time average

  		self.box1ID=self.vints["Box1Id"]
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
		
		### Defining histos
		# Event energy histogram

		self.Histos()

    ### Counters:
		self.numInputEvents = 0
		self.numOutputEvents = 0
		
		return


	############################################################
	def execute(self, event=""):

		self.m.log(2, 'Execute()')		
	
		self.numInputEvents += 1   

		#Fiducial cut. There must be one photon in the fiducial volue
		#(defined by parameters) in each box. return the interaction vertex
		#of the photon

		fiducial,self.vertexBox1,self.vertexBox2 = self.Fiducial(event)
		self.hman.fill(self.Fiducial_histo_name,fiducial)

		if fiducial != 2:
			return False

		#DT for a perfect detector

		dtTrue = self.TrueDT()
		self.m.log(2, " dtTrue =%7.2f "%(dtTrue/ps))

		self.hman.fill(self.TrueDT_histo_name,dtTrue/ps)

		
		#Compute a TimeMap including the time-ordered sequence of the
		#first PE of each SiPM hit

		self.timeMap = self.ComputeTimeMap(event)
		nhitsBox1, nhitsBox2 = self.timeMap.NumberOfHits()

		self.m.log(2, " nhits box1 = %s nhits box2 = %s"%(
			nhitsBox1,nhitsBox2))

		self.hman.fill(self.nhitsBox1_histo_name,nhitsBox1)
		self.hman.fill(self.nhitsBox2_histo_name,nhitsBox2)

		if nhitsBox1 < self.NPE or nhitsBox2 < self.NPE:
			return False
		
		self.m.log(4, " Time maps = %s"%(self.timeMap))

		#Compute DT using the first pe

		dt1 = self.DTFirstPe()

		self.m.log(2,"dt =%7.2f ps, "%(dt1/ps))
		self.hman.fill(self.DT_histo_name,dt1/ps)

		#Compute DT using up to NPE pe

		DT = self.dtNPE(dt1)

		self.m.log(2,"DT for the first %d pe"%(self.NPE))

		for i in xrange(len(DT)-1):
			DT_histo_desc = "DT"+str(i+1)
			DT_histo_name = self.alabel(DT_histo_desc)
			d = DT[i]

			self.m.log(2,"dt =%7.2f ps, "%(d/ps))
			self.hman.fill(DT_histo_name,d/ps)

		
		if self.debug == 1:
			wait()	

		self.numOutputEvents += 1
		return True

############################################################
	def finalize(self):

		self.m.log(1, 'Finalize()')

		self.m.log(1, 'Input  Events: ', self.numInputEvents)
		self.m.log(1, 'Output Events: ', self.numOutputEvents)
		

		self.logman["USER"].ints[self.alabel("InputEvents")] = self.numInputEvents
		self.logman["USER"].ints[self.alabel("OutputEvents")] = self.numOutputEvents
		

		return

############################################################
	def TrueDT(self):
		"""
		Compute the time difference between the 2 photons
		usin true information.
		"""
		t1 = self.vertexBox1.t 
		tf1 = distance((0,0,0),self.vertexBox1.XYZ())/c_light
		dt1 = t1 - tf1
		t2 = self.vertexBox2.t 
		tf2 = distance((0,0,0),self.vertexBox2.XYZ())/c_light
		dt2 = t2-tf2
		dt12 = dt1 -dt2

		self.m.log(3, " true time info")
		self.m.log(3, " t1 =%7.2f ps tf1 = %7.2f ps dt1 =%7.2f ps"%
			(t1/ps,tf1/ps,dt1/ps))
		self.m.log(3, " t2 =%7.2f ps tf2 = %7.2f ps dt2 =%7.2f ps"%
			(t2/ps,tf2/ps,dt2/ps))
		self.m.log(3, " dt12 =%7.2f "%(dt12/ps))
		
		return dt12

###########################################################
	def dtNPE(self, dt1):
		"""
		Computes the DT averaging the first NPE pes.
		"""

		self.m.log(3,
			"dtNPE: d1 =%7.2f ps "%(dt1/ps))

		dt = dt1
		DT=[]
		DT.append(dt)
		for i in xrange(self.NPE):
			dn = self.DNthPe(i+1)
			dt += dn
			self.m.log(3,
			"dtNPE: d%s =%7.2f ps "%(i+1,dn/ps))
			DT.append(dt)

		DTA=[]
		num = 1.
		for d in DT:
			d = d/num
			DTA.append(d)
			num+=1

		return DTA
###########################################################
	def DTFirstPe(self):
		"""
		Compute DT using the first PE
		"""

		siPMHit1 = self.timeMap.timeHit(boxNumber=1,index=0,jitter=self.TJ)
		siPMHit2 = self.timeMap.timeHit(boxNumber=2,index=0, jitter=self.TJ)

		self.m.log(3," DTFirstPe: Hit 1 =  %s"%(siPMHit1))
		self.m.log(3," DTFirstPe: Hit 2 =  %s"%(siPMHit2))
				
		self.hman.fill(self.T0Box1_histo_name,siPMHit1.T()/ps)
		self.hman.fill(self.T0Box2_histo_name,siPMHit2.T()/ps)

		dbox1 = distance(siPMHit1.XYZ(),self.vertexBox1.XYZ())
		tpath1 = dbox1*self.lxe.RefractionIndexUV()/c_light

		dbox2 = distance(siPMHit2.XYZ(),self.vertexBox2.XYZ())
		tpath2 = dbox2*self.lxe.RefractionIndexUV()/c_light

		self.m.log(3,
			"DTFirstPe: dbox1 =%7.2f mm,tpath1= %7.2f ps, time1 -tpath1 = %7.2f ps "%(
			dbox1/mm,tpath1/ps,(siPMHit1.T() - tpath1)/ps))
		self.m.log(3,
			"DTFirstPe: dbox2 =%7.2f mm,tpath2= %7.2f ps, time2 -tpath2 = %7.2f ps "%(
			dbox2/mm,tpath2/ps,(siPMHit2.T() - tpath2)/ps))

		self.hman.fill(self.DBox1_histo_name,dbox1/mm)
		self.hman.fill(self.TBox1_histo_name,tpath1/ps)
		self.hman.fill(self.Time1MinusTBox1_histo_name,(siPMHit1.T() - tpath1)/ps)

		self.hman.fill(self.DBox2_histo_name,dbox2/mm)
		self.hman.fill(self.TBox2_histo_name,tpath2/ps)
		self.hman.fill(self.Time2MinusTBox2_histo_name,(siPMHit2.T() - tpath2)/ps)

		tg1 = siPMHit1.T() - tpath1
		tg2 = siPMHit2.T() - tpath2
		tf1 = distance((0,0,0),self.vertexBox1.XYZ())/c_light
		tf2 = distance((0,0,0),self.vertexBox2.XYZ())/c_light

		dt1 = abs(tg1 - tf1)
		dt2 = abs(tg2 - tf2)
		dt12 = dt1 - dt2

		self.m.log(3, " +++DTFirstPe+++++")
		self.m.log(3, " tg1 =%7.2f ps tf1 = %7.2f ps dt1 =%7.2f ps"%
			(tg1/ps,tf1/ps,dt1/ps))
		self.m.log(3, " tg2 =%7.2f ps tf2 = %7.2f ps dt2 =%7.2f ps"%
			(tg2/ps,tf2/ps,dt2/ps))
		self.m.log(3, " dt12 =%7.2f "%(dt12/ps))

		return dt12

###########################################################
	def DNthPe(self,peIndex):
		"""
		Compute DT using the PE with peIndex (=0 is first PE)
		"""
		
		siPMHit1 = self.timeMap.timeHit(boxNumber=1,index=peIndex,jitter=self.TJ)
		siPMHit2 = self.timeMap.timeHit(boxNumber=2,index=peIndex, jitter=self.TJ)

		self.m.log(3," DNthPe: Hit 1 =  %s"%(siPMHit1))
		self.m.log(3," DNthPe: Hit 2 =  %s"%(siPMHit2))
				
		dbox1 = distance(siPMHit1.XYZ(),self.vertexBox1.XYZ())
		tpath1 = dbox1*self.lxe.RefractionIndexUV()/c_light

		dbox2 = distance(siPMHit2.XYZ(),self.vertexBox2.XYZ())
		tpath2 = dbox2*self.lxe.RefractionIndexUV()/c_light

		self.m.log(3,
			"DNthPe: dbox1 =%7.2f mm,tpath1= %7.2f ps, time1 -tpath1 = %7.2f ps "%(
			dbox1/mm,tpath1/ps,(siPMHit1.T() - tpath1)/ps))
		self.m.log(3,
			"DNthPe: dbox2 =%7.2f mm,tpath2= %7.2f ps, time2 -tpath2 = %7.2f ps "%(
			dbox2/mm,tpath2/ps,(siPMHit2.T() - tpath2)/ps))


		tg1 = siPMHit1.T() - tpath1
		tg2 = siPMHit2.T() - tpath2
		tf1 = distance((0,0,0),self.vertexBox1.XYZ())/c_light
		tf2 = distance((0,0,0),self.vertexBox2.XYZ())/c_light

		dt1 = abs(tg1 - tf1)
		dt2 = abs(tg2 - tf2)
		dt12 = dt1 - dt2

		self.m.log(3, " +++DNthPe+++++")
		self.m.log(3, " tg1 =%7.2f ps tf1 = %7.2f ps dt1 =%7.2f ps"%
			(tg1/ps,tf1/ps,dt1/ps))
		self.m.log(3, " tg2 =%7.2f ps tf2 = %7.2f ps dt2 =%7.2f ps"%
			(tg2/ps,tf2/ps,dt2/ps))
		self.m.log(3, " dt12 =%7.2f "%(dt12/ps))

		return dt12

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
				fid+=1
			else:
				self.m.log(2,'gamma not found in box1 or in box2')
				self.m.log(2, ' x0 =%7.2f mm, y0 =%7.2f mm, z0 =%7.2f mm '%(
			x0/mm,y0/mm,z0/mm))
				self.m.log(2, ' xf =%7.2f mm, yf =%7.2f mm, zf =%7.2f mm '%(
			x/mm,y/mm,z/mm))
				self.m.log(2, "fBox1 --", self.fbox1)
				self.m.log(2, "fBox2 ---", self.fbox2)
				
				
		return (fid,vertexBox1,vertexBox2)
############################################################
	def ComputeTimeMap(self,event):
		"""
		Compute time maps
		"""

		sensorhits =  event.GetMCSensHits()
		self.m.log(3, " event has %d sensor hits = "%(len(sensorhits)))

		timeMapBox1={}
		timeMapBox2={}

		for hit in sensorhits:
			hid = hit.GetSensorID()
			xh = hit.GetPosition().x()
			yh = hit.GetPosition().y()
			zh = hit.GetPosition().z()
			Ah = hit.GetAmplitude()
			self.m.log(4, " hit, ID = ", hid)
			self.m.log(4, ' xh =%7.2f mm, yh =%7.2f mm, zh =%7.2f mm Q =%d pes '%(
			xh/mm,yh/mm,zh/mm,Ah))

			self.m.log(4, " waveform for hit ID = %d"%(hit.GetSensorID()))
			waveform = hit.GetWaveform().GetData()
			
			for timeBins in waveform:
				# get the arrival time of the pe
				tbin = timeBins.first
				time = tbin*5*ps
				A = timeBins.second

				self.m.log(4, " tbin = %d time = %7.2f ps A = %7.2f pes"%(
					tbin,time/ps,A))
				#A is the number of pes in the bin. It is equal to 1
				# most of the time (because the time bining id very thin)
				#but it can be a higher number. A is computed assuming QE=1
				# One needs to correct for the QE of the device. 

				npes = 0
				for p in xrange(0,A):
					if rnd.uniform(0.,1.) <= self.QE:
						npes+=1

				self.m.log(4, " npes = %d "%(
					npes))
				#is there a real pe in this time bin? otherwise loop
				if npes < 1:
					continue
				else:
					sipm=[]
					sipm.append(xh)
					sipm.append(yh)
					sipm.append(zh)
					sipm.append(Ah)
					sipm.append(time)
					self.m.log(4, " ID = %d time= %7.2f "%(
						hid, time/ps))

					if self.boxId(hid) == 1:
						timeMapBox1[hid] = sipm

					elif self.boxId(hid) == 2:
						timeMapBox2[hid] = sipm 
					else:
						print "error: hit index not in range, hid=%d"%(hid)
						sys.exit(0)

					break;

		# sort the maps according to the time stamp of pes
		TimeMapBox1 = sorted(timeMapBox1.items(), key=sortHits)
		TimeMapBox2 = sorted(timeMapBox2.items(), key=sortHits)

		self.m.log(4, " sorted time map box1 = ",TimeMapBox1)
		self.m.log(4, " sorted time map box2 = ",TimeMapBox2)

		timeMap = TimeMap(TimeMapBox1,TimeMapBox2)

		return timeMap
		

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
			print i
			lbl = blb+str(i)
			print blb
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
	def PrintTimes(self,TimeMapBox):
		
		s=''
		for hit in TimeMapBox:
			hid = hit[0]
			values = hit[1]
			time = values[4]/ps
			s+="ID = %s time =%7.2f ps \n"%(hid,time)

		return s

	############################################################		
	def Histos(self):
		"""
		book the histograms for the algo 
		"""
		TrueDT_histo_desc = "TrueDT"
		self.TrueDT_histo_name = self.alabel(TrueDT_histo_desc)
		self.hman.h1(self.TrueDT_histo_name, TrueDT_histo_desc, 
			20, -5, 5)
		self.hman.fetch(
			self.TrueDT_histo_name).GetXaxis().SetTitle(
			"True DT box1 - box2 ps")

		Fiducial_histo_desc = "Fiducial"
		self.Fiducial_histo_name = self.alabel(Fiducial_histo_desc)
		self.hman.h1(self.Fiducial_histo_name, Fiducial_histo_desc, 
			10, 0, 5)
		self.hman.fetch(
			self.Fiducial_histo_name).GetXaxis().SetTitle(
			"Fiducial interactions")

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

		T0Box1_histo_desc = "T0Box1"
		self.T0Box1_histo_name = self.alabel(T0Box1_histo_desc)
		self.hman.h1(self.T0Box1_histo_name, T0Box1_histo_desc, 
			50, 300, 800)
		self.hman.fetch(
			self.T0Box1_histo_name).GetXaxis().SetTitle(
			"t0 box1 (ps)")

		DBox1_histo_desc = "DBox1"
		self.DBox1_histo_name = self.alabel(DBox1_histo_desc)
		self.hman.h1(self.DBox1_histo_name, DBox1_histo_desc, 
			50, 0, 50)
		self.hman.fetch(
			self.DBox1_histo_name).GetXaxis().SetTitle(
			"distance T0 hit to vertex box 1 (mm)")

		TBox1_histo_desc = "TBox1"
		self.TBox1_histo_name = self.alabel(TBox1_histo_desc)
		self.hman.h1(self.TBox1_histo_name, TBox1_histo_desc, 
			50, 0, 500)
		self.hman.fetch(
			self.TBox1_histo_name).GetXaxis().SetTitle(
			"time T0 hit to vertex box 1 (ps)")

		DBox2_histo_desc = "DBox2"
		self.DBox2_histo_name = self.alabel(DBox2_histo_desc)
		self.hman.h1(self.DBox2_histo_name, DBox2_histo_desc, 
			50, 0, 50)
		self.hman.fetch(
			self.DBox2_histo_name).GetXaxis().SetTitle(
			"distance T0 hit to vertex box 1 (mm)")

		TBox2_histo_desc = "TBox2"
		self.TBox2_histo_name = self.alabel(TBox2_histo_desc)
		self.hman.h1(self.TBox2_histo_name, TBox2_histo_desc, 
			50, 0, 500)
		self.hman.fetch(
			self.TBox2_histo_name).GetXaxis().SetTitle(
			"time T0 hit to vertex box 2 (ps)")

		Time1MinusTBox1_histo_desc = "Time1MinusTBox1"
		self.Time1MinusTBox1_histo_name = self.alabel(Time1MinusTBox1_histo_desc)
		self.hman.h1(self.Time1MinusTBox1_histo_name, Time1MinusTBox1_histo_desc, 
			50, -500, 500)
		self.hman.fetch(
			self.Time1MinusTBox1_histo_name).GetXaxis().SetTitle(
			"T0 minus tvertex box1 (ps)")

		Time2MinusTBox2_histo_desc = "Time2MinusTBox2"
		self.Time2MinusTBox2_histo_name = self.alabel(Time2MinusTBox2_histo_desc)
		self.hman.h1(self.Time2MinusTBox2_histo_name, Time2MinusTBox2_histo_desc, 
			50, -500, 500)
		self.hman.fetch(
			self.Time2MinusTBox2_histo_name).GetXaxis().SetTitle(
			"T0 minus tvertex box2 (ps)")

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

		T0Box2_histo_desc = "T0Box2"
		self.T0Box2_histo_name = self.alabel(T0Box2_histo_desc)
		self.hman.h1(self.T0Box2_histo_name, T0Box2_histo_desc, 
			50,300, 800)

		self.hman.fetch(
			self.T0Box2_histo_name).GetXaxis().SetTitle(
			"t0 box2 (ps)")


		DT_histo_desc = "DT"
		self.DT_histo_name = self.alabel(DT_histo_desc)
		self.hman.h1(self.DT_histo_name, DT_histo_desc, 
			50, -200, 200)
		self.hman.fetch(
			self.DT_histo_name).GetXaxis().SetTitle(
			"DT (ps)")

		for i in xrange(self.NPE):
			DT_histo_desc = "DT"+str(i+1)
			lbl ="DT computed with %s PES (ps)"%(i+1)
			DT_histo_name = self.alabel(DT_histo_desc)
			self.hman.h1(DT_histo_name, DT_histo_desc, 
			50, -200, 200)
			self.hman.fetch(DT_histo_name).GetXaxis().SetTitle(lbl)




	

	
	