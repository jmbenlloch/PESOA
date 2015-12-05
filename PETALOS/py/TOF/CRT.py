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

############################################################
def sortHits(hit):
	hitVal = hit[1]
	return hitVal[4]

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
  		self.QE = self.doubles["QE"]
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

		# self.lxe = LXe() #lxe properties
		# print self.lxe

		if self.debug == 1:
			wait()
   

	############################################################		
	def initialize(self):

		self.m.log(1, 'Initialize()')
		
		### Defining histos
		# Event energy histogram

		self.XYZBox1_histo_desc = "XYZBox1"
		self.XYZBox1_histo_name = self.alabel(self.XYZBox1_histo_desc)
		self.hman.h3(self.XYZBox1_histo_name, self.XYZBox1_histo_desc, 
			10, self.box1.xmin, self.box1.xmax,
           	10, self.box1.ymin, self.box1.ymax,
           	10, self.box1.zmin, self.box1.zmax)
		self.hman.fetch(
			self.XYZBox1_histo_name).GetXaxis().SetTitle(
			"XYZ interaction box 1 (mm)")

		self.XYBox1_histo_desc = "XYBox1"
		self.XYBox1_histo_name = self.alabel(self.XYBox1_histo_desc)
		self.hman.h2(self.XYBox1_histo_name, self.XYBox1_histo_desc, 
			25, self.box1.xmin, self.box1.xmax,
           	25, self.box1.ymin, self.box1.ymax)  	
		self.hman.fetch(
			self.XYBox1_histo_name).GetXaxis().SetTitle(
			"XY interaction box 1 (mm)")

		self.ZBox1_histo_desc = "ZBox1"
		self.ZBox1_histo_name = self.alabel(self.ZBox1_histo_desc)
		self.hman.h1(self.ZBox1_histo_name, self.ZBox1_histo_desc, 
			50, self.box1.zmin, self.box1.zmax)
		self.hman.fetch(
			self.ZBox1_histo_name).GetXaxis().SetTitle(
			"Z interaction box 1 (mm)")

		self.T0Box1_histo_desc = "T0Box1"
		self.T0Box1_histo_name = self.alabel(self.T0Box1_histo_desc)
		self.hman.h1(self.T0Box1_histo_name, self.T0Box1_histo_desc, 
			50, 0, 500)
		self.hman.fetch(
			self.T0Box1_histo_name).GetXaxis().SetTitle(
			"t0 box1 (ps)")

		self.DBox1_histo_desc = "DBox1"
		self.DBox1_histo_name = self.alabel(self.DBox1_histo_desc)
		self.hman.h1(self.DBox1_histo_name, self.DBox1_histo_desc, 
			50, 0, 50)
		self.hman.fetch(
			self.DBox1_histo_name).GetXaxis().SetTitle(
			"distance T0 hit to vertex box 1 (mm)")

		self.TBox1_histo_desc = "TBox1"
		self.TBox1_histo_name = self.alabel(self.TBox1_histo_desc)
		self.hman.h1(self.TBox1_histo_name, self.TBox1_histo_desc, 
			50, 0, 500)
		self.hman.fetch(
			self.TBox1_histo_name).GetXaxis().SetTitle(
			"time T0 hit to vertex box 1 (ps)")

		self.DBox2_histo_desc = "DBox2"
		self.DBox2_histo_name = self.alabel(self.DBox2_histo_desc)
		self.hman.h1(self.DBox2_histo_name, self.DBox2_histo_desc, 
			50, 0, 50)
		self.hman.fetch(
			self.DBox2_histo_name).GetXaxis().SetTitle(
			"distance T0 hit to vertex box 1 (mm)")

		self.TBox2_histo_desc = "TBox2"
		self.TBox2_histo_name = self.alabel(self.TBox2_histo_desc)
		self.hman.h1(self.TBox2_histo_name, self.TBox2_histo_desc, 
			50, 0, 500)
		self.hman.fetch(
			self.TBox2_histo_name).GetXaxis().SetTitle(
			"time T0 hit to vertex box 2 (ps)")

		self.Time1MinusTBox1_histo_desc = "Time1MinusTBox1"
		self.Time1MinusTBox1_histo_name = self.alabel(self.Time1MinusTBox1_histo_desc)
		self.hman.h1(self.Time1MinusTBox1_histo_name, self.Time1MinusTBox1_histo_desc, 
			50, -200, 200)
		self.hman.fetch(
			self.Time1MinusTBox1_histo_name).GetXaxis().SetTitle(
			"T0 minus tvertex box1 (ps)")

		self.Time2MinusTBox2_histo_desc = "Time2MinusTBox2"
		self.Time2MinusTBox2_histo_name = self.alabel(self.Time2MinusTBox2_histo_desc)
		self.hman.h1(self.Time2MinusTBox2_histo_name, self.Time2MinusTBox2_histo_desc, 
			50, -200, 200)
		self.hman.fetch(
			self.Time2MinusTBox2_histo_name).GetXaxis().SetTitle(
			"T0 minus tvertex box2 (ps)")

			

		self.XYZBox2_histo_desc = "XYZBox2"
		self.XYZBox2_histo_name = self.alabel(self.XYZBox2_histo_desc)
		self.hman.h3(self.XYZBox2_histo_name, self.XYZBox2_histo_desc, 
			10, self.box2.xmin, self.box2.xmax,
           	10, self.box2.ymin, self.box2.ymax,
           	10, self.box2.zmin, self.box2.zmax)
		self.hman.fetch(
			self.XYZBox2_histo_name).GetXaxis().SetTitle(
			"XYZ interaction box 2 (mm)")

		self.XYBox2_histo_desc = "XYBox2"
		self.XYBox2_histo_name = self.alabel(self.XYBox2_histo_desc)
		self.hman.h2(self.XYBox2_histo_name, self.XYBox2_histo_desc, 
			25, self.box2.xmin, self.box2.xmax,
           	25, self.box2.ymin, self.box2.ymax)
           	
		self.hman.fetch(
			self.XYBox2_histo_name).GetXaxis().SetTitle(
			"XY interaction box 2 (mm)")

		self.ZBox2_histo_desc = "ZBox2"
		self.ZBox2_histo_name = self.alabel(self.ZBox2_histo_desc)
		self.hman.h1(self.ZBox2_histo_name, self.ZBox2_histo_desc, 
			50, self.box1.zmin, self.box1.zmax)
   
		self.hman.fetch(
			self.ZBox2_histo_name).GetXaxis().SetTitle(
			"Z interaction box 2 (mm)")

		self.T0Box2_histo_desc = "T0Box2"
		self.T0Box2_histo_name = self.alabel(self.T0Box2_histo_desc)
		self.hman.h1(self.T0Box2_histo_name, self.T0Box2_histo_desc, 
			50, 0, 500)

		self.hman.fetch(
			self.T0Box2_histo_name).GetXaxis().SetTitle(
			"t0 box2 (ps)")

		self.T0Box12_histo_desc = "T0Box12"
		self.T0Box12_histo_name = self.alabel(self.T0Box12_histo_desc)
		self.hman.h2(self.T0Box12_histo_name, self.T0Box12_histo_desc, 
			20, 0, 200,
           	20, 0, 200)

		self.hman.fetch(
			self.T0Box12_histo_name).GetXaxis().SetTitle(
			"t0 box1 vs box2 (ps)")

		self.T0Box1MinusT0Box2_histo_desc = "T0Box1MinusT0Box2"
		self.T0Box1MinusT0Box2_histo_name = self.alabel(self.T0Box1MinusT0Box2_histo_desc)
		self.hman.h1(self.T0Box1MinusT0Box2_histo_name, self.T0Box1MinusT0Box2_histo_desc, 
			50, -250, 250)
		self.hman.fetch(
			self.T0Box1MinusT0Box2_histo_name).GetXaxis().SetTitle(
			"t0 box1 minus t0 box2 (ps)")

		self.DT_histo_desc = "DT"
		self.DT_histo_name = self.alabel(self.DT_histo_desc)
		self.hman.h1(self.DT_histo_name, self.DT_histo_desc, 
			50, -200, 200)
		self.hman.fetch(
			self.DT_histo_name).GetXaxis().SetTitle(
			"DT (ps)")
    ### Counters:
		self.numInputEvents = 0
		self.numOutputEvents = 0
		
		return



	############################################################
	def execute(self, event=""):

		self.m.log(2, 'Execute()')		
	
		self.numInputEvents += 1   
		primaryParticles = PrimaryParticles(event)
		self.m.log(2, ' number of primary Particles =%d '%(len(primaryParticles)))
		
		self.vertexBox1=Point3D()
		self.vertexBox2=Point3D()
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

				self.vertexBox1.x = x
				self.vertexBox1.y = y
				self.vertexBox1.z = z
			
			elif self.fbox2.Active((x,y,z)) == True:
				self.m.log(2,'gamma found in box2')
				
				self.hman.fill(self.XYZBox2_histo_name, 
					x/mm,y/mm,z/mm)
				self.hman.fill(self.XYBox2_histo_name, 
					x/mm,y/mm)
				self.hman.fill(self.ZBox2_histo_name, 
					z/mm)

				self.vertexBox2.x = x
				self.vertexBox2.y = y
				self.vertexBox2.z = z
			else:
				self.m.log(2,'gamma not found in box1 or in box2')
				self.m.log(2, ' x0 =%7.2f mm, y0 =%7.2f mm, z0 =%7.2f mm '%(
			x0/mm,y0/mm,z0/mm))
				self.m.log(2, ' xf =%7.2f mm, yf =%7.2f mm, zf =%7.2f mm '%(
			x/mm,y/mm,z/mm))
				self.m.log(2, "fBox1 --", self.fbox1)
				self.m.log(2, "fBox2 ---", self.fbox2)
				
				return False

			
		self.numOutputEvents += 1

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
			self.m.log(3, " hit, ID = ", hid)
			self.m.log(3, ' xh =%7.2f mm, yh =%7.2f mm, zh =%7.2f mm Q =%d pes  '%(
			xh/mm,yh/mm,zh/mm,Ah))

			self.m.log(4, " waveform for hit ID = %d"%(hit.GetSensorID()))
			waveform = hit.GetWaveform().GetData()
			
			for timeBins in waveform:
				# get the arrival time of the first pe
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

		TimeMapBox1 = sorted(timeMapBox1.items(), key=sortHits)
		TimeMapBox2 = sorted(timeMapBox2.items(), key=sortHits)

		self.m.log(4, " sorted time map box1 = ",TimeMapBox1)
		self.m.log(4, " sorted time map box2 = ",TimeMapBox2)

		timeMap = TimeMap(TimeMapBox1,TimeMapBox2)
		
		self.m.log(3, " Time maps = %s"%(timeMap))

		x1,y1,z1,A1,time1 =timeMap.timeHit(boxNumber=1,index=0)
		x2,y2,z2,A2,time2 =timeMap.timeHit(boxNumber=2,index=0)

		self.m.log(2,
			"timeHit 1:xh =%7.2f mm,yh =%7.2f mm,zh =%7.2f mm, A = %7.2f pes t= %7.2f ps"%(
			x1/mm,y1/mm,z1/mm,A1,time1/ps))
		self.m.log(2,
			"timeHit 2:xh =%7.2f mm,yh =%7.2f mm,zh =%7.2f mm, A = %7.2f pes t= %7.2f ps"%(
			x2/mm,y2/mm,z2/mm,A2,time2/ps))

		self.hman.fill(self.T0Box1_histo_name,time1/ps)
		self.hman.fill(self.T0Box2_histo_name,time2/ps)
		self.hman.fill(self.T0Box12_histo_name,time1/ps,time2/ps)
		self.hman.fill(self.T0Box1MinusT0Box2_histo_name,(time1-time2)/ps) 
			

		dbox1 = distance((x1,y1,z1),self.vertexBox1.XYZ())
		tpath1 = dbox1/c_light

		dbox2 = distance((x2,y2,z2),self.vertexBox2.XYZ())
		tpath2 = dbox2/c_light

		self.m.log(2,
			"dbox1 =%7.2f mm,tpath1= %7.2f ps, time1 -tpath1 = %7.2f ps "%(
			dbox1/mm,tpath1/ps,(time1 - tpath1)/ps))
		self.m.log(2,
			"dbox2 =%7.2f mm,tpath2= %7.2f ps, time2 -tpath2 = %7.2f ps "%(
			dbox2/mm,tpath2/ps,(time2 - tpath2)/ps))

		self.hman.fill(self.DBox1_histo_name,dbox1/mm)
		self.hman.fill(self.TBox1_histo_name,tpath1/ps)
		self.hman.fill(self.Time1MinusTBox1_histo_name,(time1 - tpath1)/ps)

		self.hman.fill(self.DBox2_histo_name,dbox2/mm)
		self.hman.fill(self.TBox2_histo_name,tpath2/ps)
		self.hman.fill(self.Time2MinusTBox2_histo_name,(time2 - tpath2)/ps)

		dt1 = time1 - tpath1
		dt2 = time2 - tpath2
		dt = dt1 - dt2

		self.m.log(2,"dt =%7.2f ps, "%(dt/ps))
		self.hman.fill(self.DT_histo_name,dt/ps)
		

		if self.debug == 1:
			wait()	
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
			self.m.log(1, "loading parameter", lbl)
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



	

	
	