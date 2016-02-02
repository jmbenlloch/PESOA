from Centella.AAlgo import AAlgo
from Util import *
from Geometry import *
from TOF import *
from Scintillator import *
import random as rnd


"""
This algorithm computes the Time map. 
"""


class ATimeMap(AAlgo):

	############################################################
	def __init__(self, param=False, level=1, label="", **kargs):

		"""
		ATimeMap Algorithm
		"""
		#self.m.log(1, 'Constructor()')

		### GENERAL STUFF
		self.name = 'ATimeMap'
		
		AAlgo.__init__(self, param, level, self.name, 0, label, kargs)

    ### PARAMETERS
		
		self.debug = self.ints["Debug"]  #used to stop program at key break points

  		self.SCINT = self.strings["SCINTILLATOR"]
  		self.QE = self.doubles["QE"]  #quantum efficiency
  		self.DTMAX = self.doubles["DTMAX"]*ps  #max diff wrt first pes
  		self.SPTR = self.doubles["SPTR"]*ps  #single photon time resolution
  		self.ASIC= self.doubles["ASIC"]*ps  #ASIC contribution
  		self.INTER = self.strings["INTER"]
 	
  		self.NSIPM = self.ints["NSIPM"]   #number of SiPMs per box

  		#time jitter of SiPM + ASIC
  		self.TJ = sqrt(self.SPTR**2+self.ASIC**2)
  		self.box1ID=self.vints["Box1Id"]  #ids of SiPMs in box1
  		self.box2ID=self.vints["Box2Id"] 

  		self.HBOX = self.ints["HBOX"] # if 1 only box1 histos if 2 box1 and 2 if 0 no box histos
		self.m.log(1, "IDs Box1 --", self.box1ID)
		self.m.log(1, "IDs Box2 --", self.box2ID)
		self.m.log(1, "Scintillator --", self.SCINT)

		self.m.log(1, "QE = %7.2f Time Jitter =%7.2f ps  --"%(self.QE, self.TJ/ps))

		if self.SCINT == "LXE":
			self.scint = LXe() #lxe properties
		elif self.SCINT == "LYSO":
			self.scint = LYSO()
		else:
			print "scintillator not yet implemented"
			sys.exit()

		print self.scint

		
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

	#Time Map instance
		self.timeMap = TimeMap(numberOfBoxes = 2, dtmax=self.DTMAX)
		
		return


	############################################################
	def execute(self, event=""):

		self.m.log(2, 'Execute()')		
	
		self.numInputEvents += 1   

		
		#Compute a TimeMap including the time-ordered sequence of the
		#PES of each SiPM hit that arrive within DTMAX of the first PE

		
		
		self.timeMap.SetInteractionVertices(self.logman["USER"].gparam["vertexBox"])

		self.ComputeTimeMap(event)
		nhitsBox1 = self.timeMap.NumberOfSiPmHits(box=1)
		nhitsBox2 = self.timeMap.NumberOfSiPmHits(box=2)

		self.m.log(2, " nhits box1 = %s nhits box2 = %s"%(
			nhitsBox1,nhitsBox2))

		if self.HBOX ==1 :
			self.hman.fill(self.nhitsBox1_histo_name,nhitsBox1)
		if self.HBOX ==2 :
			self.hman.fill(self.nhitsBox2_histo_name,nhitsBox2)

		if nhitsBox1 < self.NSIPM or nhitsBox2 < self.NSIPM:
			return False
		
		self.m.log(3, " Time maps = %s"%(self.timeMap))
		if self.debug == 1:
			wait()	
		
		self.numOutputEvents += 1
		self.logman["USER"].gparam["TimeMap"]=self.timeMap
		return True

############################################################
	def finalize(self):

		self.m.log(1, 'Finalize()')

		self.m.log(1, 'Input  Events: ', self.numInputEvents)
		self.m.log(1, 'Output Events: ', self.numOutputEvents)

		return



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

		vertexBox1 =self.timeMap.InteractionVertex(box=1)
		vertexBox2 =self.timeMap.InteractionVertex(box=2)

		ng1 = 0
		ng2 = 0

#                for w in sensorhits[0].GetWaveform().GetData():
#                    print "%s \t %s" % (w.first,w.second)

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

			sipmhit = SiPMHit(hid,xh,yh,zh,Ah,time0,
				self.QE,self.SPTR,self.ASIC,self.DTMAX)

			self.m.log(5, " ++++++hit = %s"%(sipmhit))
			
			if self.debug == 2:
				wait()
			
			if self.boxId(hid) == 1:
				ng1+=Ah
			elif self.boxId(hid) == 2:
				ng2+=Ah
			else:
				print "error: hit index not in range, hid=%d"%(hid)
				sys.exit(0)
			
			dbox1 = distance(sipmhit.XYZ(),vertexBox1.XYZ())
			tpath1 = dbox1*self.scint.RefractionIndex()/c_light
			dbox2 = distance(sipmhit.XYZ(),vertexBox2.XYZ())
			tpath2 = dbox2*self.scint.RefractionIndex()/c_light

			# keep all the pes within DTMAX (~300 ps) of first pe
			np=0
			Q=0
			timeFirstPe = 1e+9*ns

			for timeBins in waveform:
				# get the arrival time of the pe
				tbin = timeBins.first
				time = tbin*5*ps
				A = timeBins.second
				np+=1

				self.m.log(6, "pe number = %d tbin = %d time = %7.2f ps A = %7.2f pes"%(
					np, tbin,time/ps,A))
				self.m.log(6, " DT wrt 1st = %7.2f ps"%((time - time0)/ps))

				if self.boxId(hid) == 1:
					tpes  = time - vertexBox1.t -tpath1
					self.m.log(6, "BOX1: pe subtracted time = %7.2f ps"%(tpes/ps))

					if self.HBOX ==1 :
						self.hman.fill(self.TpesBOX1_histo_name,tpes/ps)
						self.hman.fill(self.DTpesBOX1_histo_name,(time - time0)/ps)
					if len(sipmhit.W)==0:
						self.m.log(5, "BOX1: 1st PE subtracted time = %7.2f ps"%(tpes/ps))
						if self.HBOX ==1 :
							self.hman.fill(self.T0pesBOX1_histo_name,tpes/ps)
				else:
					tpes  = time - vertexBox2.t -tpath2
					self.m.log(6, "BOX2: pe subtracted time = %7.2f ps"%(tpes/ps))
					if self.HBOX ==2 :
						self.hman.fill(self.TpesBOX2_histo_name,tpes/ps)
						self.hman.fill(self.DTpesBOX2_histo_name,(time - time0)/ps)
					if len(sipmhit.W)==0:
						self.m.log(5, "BOX2: 1st PE subtracted time = %7.2f ps"%(tpes/ps))
						if self.HBOX ==2 :
							self.hman.fill(self.T0pesBOX2_histo_name,tpes/ps)
					
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
				else: 
					#valid pes, 

					if self.boxId(hid) == 1:
						tpes  = time - vertexBox1.t -tpath1
						self.m.log(6, "BOX1: pe subtracted time after QE = %7.2f ps"%(tpes/ps))
						if self.HBOX ==1 :
							if self.QE <1 :
								self.hman.fill(self.TpesQeBOX1_histo_name,tpes/ps)
						if len(sipmhit.W)==0:
							self.m.log(5, "BOX1: 1st pe subtracted time after QE = %7.2f ps"%(tpes/ps))
							if self.HBOX ==1 :
								if self.QE <1 :
									self.hman.fill(self.T0pesQeBOX1_histo_name,tpes/ps)
					else:
						tpes  = time - vertexBox2.t -tpath2
						self.m.log(6, "BOX2: pe subtracted time after QE = %7.2f ps"%(tpes/ps))
						if self.HBOX ==2 :
							if self.QE <1 :
								self.hman.fill(self.TpesQeBOX2_histo_name,tpes/ps)
						if len(sipmhit.W)==0:
							self.m.log(5, "BOX2: 1st pe subtracted time after QE = %7.2f ps"%(tpes/ps))
							if self.HBOX ==2 :
								if self.QE <1 :
									self.hman.fill(self.T0pesQeBOX2_histo_name,tpes/ps)

					# smear time, to take into account the effect of ASIC and SiPM
					stime = SmearTime(time,self.SPTR,self.ASIC)

					if self.boxId(hid) == 1:
		
						tpes  = stime - vertexBox1.t -tpath1
						self.m.log(6, "BOX1: pe subtracted time after QE and Smear = %7.2f ps"%(
							tpes/ps))
						if self.HBOX ==1 :
							if self.SPTR > 0 or self.ASIC > 0 :
								self.hman.fill(self.TpesQeSmearBOX1_histo_name,tpes/ps)
						if len(sipmhit.W)==0:
							self.m.log(6, "BOX1: 1st pe subtracted time after QE and Smear = %7.2f ps"%(
							tpes/ps))
							if self.HBOX ==1 :
								if self.SPTR > 0 or self.ASIC > 0 :
									self.hman.fill(self.T0pesQeSmearBOX1_histo_name,tpes/ps)
					else:
						tpes  = stime - vertexBox2.t -tpath2
						self.m.log(6, "BOX2: pe subtracted time after QE and Smear= %7.2f ps"%(
							tpes/ps))
						if self.HBOX ==2 :
							if self.SPTR > 0 or self.ASIC > 0 :
								self.hman.fill(self.TpesQeSmearBOX2_histo_name,tpes/ps)
						if len(sipmhit.W)==0:
							self.m.log(6, "BOX2: 1st pe subtracted time after QE and Smear= %7.2f ps"%(
							tpes/ps))
							if self.HBOX ==2 :
								if self.SPTR > 0 or self.ASIC > 0 :
									self.hman.fill(self.T0pesQeSmearBOX2_histo_name,tpes/ps)

					#keep the time stamp if forst PE ortime within DTMAX

					if len(sipmhit.W)==0: #first PE
						sipmhit.W.append(stime)
						timeFirstPe = stime 

					elif abs(stime - timeFirstPe) < self.DTMAX:
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
				if self.HBOX ==1 :
					self.hman.fill(self.npesBOX1_histo_name,sipmhit.NumberOfPE())

			elif self.boxId(hid) == 2:
				timeMapBox2[hid] = sipmhit 
				if self.HBOX ==2 :
					self.hman.fill(self.npesBOX2_histo_name,sipmhit.NumberOfPE())
			else:
				print "error: hit index not in range, hid=%d"%(hid)
				sys.exit(0)
			
		if self.HBOX ==1 :
			self.hman.fill(self.NGBOX1_histo_name,ng1)
		if self.HBOX ==2 :	
			self.hman.fill(self.NGBOX2_histo_name,ng2)
		
		if ng1 ==0 or ng2 ==0:
			print "!!! ng1 = %d, ng2=%d"%(ng1,ng2)
			return False

		if len(timeMapBox1) ==0 or len(timeMapBox2) ==0:
			print "!!! len(timeMapBox1) = %d, len(timeMapBox2)=%d"%(
				len(timeMapBox1),len(timeMapBox2))
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
	def Histos(self):
		"""
		book the histograms for the algo 
		"""

		if self.HBOX ==1 :
			NGBOX1_histo_desc = "NGBOX1"
			self.NGBOX1_histo_name = self.alabel(NGBOX1_histo_desc)

			if self.INTER == "CHER":
				self.hman.h1(self.NGBOX1_histo_name, NGBOX1_histo_desc, 
				25, 0, 200)
			else:
				self.hman.h1(self.NGBOX1_histo_name, NGBOX1_histo_desc, 
				100, 0, 30000)

			self.hman.fetch(
			self.NGBOX1_histo_name).GetXaxis().SetTitle(
			"Number of gammas in box1")

			npesBOX1_histo_desc = "NPESBOX1"
			self.npesBOX1_histo_name = self.alabel(npesBOX1_histo_desc)

			if self.INTER == "CHER":
				self.hman.h1(self.npesBOX1_histo_name, npesBOX1_histo_desc, 
				25, 0, 50)
			else:
				self.hman.h1(self.npesBOX1_histo_name, npesBOX1_histo_desc, 
				100, 0, 500)
			self.hman.fetch(
			self.npesBOX1_histo_name).GetXaxis().SetTitle(
			"Number of pes/SiPM in box1")

			nhitsBox1_histo_desc = "nhitsBox1"
			self.nhitsBox1_histo_name = self.alabel(nhitsBox1_histo_desc)
			self.hman.h1(self.nhitsBox1_histo_name, nhitsBox1_histo_desc, 
			100, 0, 200)
			self.hman.fetch(
			self.nhitsBox1_histo_name).GetXaxis().SetTitle(
			"Number of sensor hits in box1")

			T0pesBOX1_histo_desc = "TimeFirstPESubtractedBox1"
			self.T0pesBOX1_histo_name = self.alabel(T0pesBOX1_histo_desc)
			self.hman.h1(self.T0pesBOX1_histo_name, T0pesBOX1_histo_desc, 
			100, 0, 1000)
			self.hman.fetch(
			self.T0pesBOX1_histo_name).GetXaxis().SetTitle(
			"Time FirstPE Subtracted Box1 (ps)")

			TpesBOX1_histo_desc = "TimePESubtractedBox1"
			self.TpesBOX1_histo_name = self.alabel(TpesBOX1_histo_desc)
			self.hman.h1(self.TpesBOX1_histo_name, TpesBOX1_histo_desc, 
			100, 0, 1000)
			self.hman.fetch(
			self.TpesBOX1_histo_name).GetXaxis().SetTitle(
			"Time Inclusive PE Subtracted Box1 (ps)")

			DTpesBOX1_histo_desc = "TimeFirstPeMinusTimePeBox1"
			self.DTpesBOX1_histo_name = self.alabel(DTpesBOX1_histo_desc)
			self.hman.h1(self.DTpesBOX1_histo_name, DTpesBOX1_histo_desc, 
			100, 0, 1000)
			self.hman.fetch(
			self.DTpesBOX1_histo_name).GetXaxis().SetTitle(
			"Time Inclusive PE minus first PE Box1 (ps)")

			if self.QE <1:

				T0pesQeBOX1_histo_desc = "TimeFirstPESubtractedAfterQEBox1"
				self.T0pesQeBOX1_histo_name = self.alabel(T0pesQeBOX1_histo_desc)
				self.hman.h1(self.T0pesQeBOX1_histo_name, T0pesQeBOX1_histo_desc, 
				100, 0, 1000)
				self.hman.fetch(
				self.T0pesQeBOX1_histo_name).GetXaxis().SetTitle(
				"Time FirstPE Subtracted Box1 after QE (ps)")

				TpesQeBOX1_histo_desc = "TimePESubtractedAfterQEBox1"
				self.TpesQeBOX1_histo_name = self.alabel(TpesQeBOX1_histo_desc)
				self.hman.h1(self.TpesQeBOX1_histo_name, TpesQeBOX1_histo_desc, 
				100, 0, 1000)
				self.hman.fetch(
				self.TpesQeBOX1_histo_name).GetXaxis().SetTitle(
				"Time Inclusive PE Subtracted after QE Box1 (ps)")

			if self.SPTR > 0 or self.ASIC > 0 :
				T0pesQeSmearBOX1_histo_desc = "TimeFirstPESubtractedAfterQESmearBox1"
				self.T0pesQeSmearBOX1_histo_name = self.alabel(T0pesQeSmearBOX1_histo_desc)
				self.hman.h1(self.T0pesQeSmearBOX1_histo_name, T0pesQeSmearBOX1_histo_desc, 
				100, 0, 1000)
				self.hman.fetch(
				self.T0pesQeSmearBOX1_histo_name).GetXaxis().SetTitle(
				"Time FirstPE Subtracted Box1 after QE and smear (ps)")

				TpesQeSmearBOX1_histo_desc = "TimePESubtractedAfterQESmearBox1"
				self.TpesQeSmearBOX1_histo_name = self.alabel(TpesQeSmearBOX1_histo_desc)
				self.hman.h1(self.TpesQeSmearBOX1_histo_name, TpesQeSmearBOX1_histo_desc, 
				100, 0, 1000)
				self.hman.fetch(
				self.TpesQeSmearBOX1_histo_name).GetXaxis().SetTitle(
				"Time Inclusive PE Subtracted after QE and smear Box1 (ps)")

			

		if self.HBOX ==2 :
			NGBOX2_histo_desc = "NGBOX2"
			self.NGBOX2_histo_name = self.alabel(NGBOX2_histo_desc)

			if self.INTER == "CHER":
				self.hman.h1(self.NGBOX2_histo_name, NGBOX2_histo_desc, 
				25, 0, 200)
			else:
				self.hman.h1(self.NGBOX2_histo_name, NGBOX2_histo_desc, 
				100, 0, 30000)

			self.hman.fetch(
			self.NGBOX2_histo_name).GetXaxis().SetTitle(
			"Number of gammas in box2")

			npesBOX2_histo_desc = "NPESBOX2"
			self.npesBOX2_histo_name = self.alabel(npesBOX2_histo_desc)

			if self.INTER == "CHER":
				self.hman.h1(self.npesBOX2_histo_name, npesBOX2_histo_desc, 
				25, 0, 50)
			else:
				self.hman.h1(self.npesBOX2_histo_name, npesBOX2_histo_desc, 
				100, 0, 500)
			self.hman.fetch(
			self.npesBOX2_histo_name).GetXaxis().SetTitle(
			"Number of pes/SiPM in box2")

			nhitsBox2_histo_desc = "nhitsBox2"
			self.nhitsBox2_histo_name = self.alabel(nhitsBox2_histo_desc)
			self.hman.h1(self.nhitsBox2_histo_name, nhitsBox2_histo_desc, 
			100, 0, 200)
			self.hman.fetch(
			self.nhitsBox2_histo_name).GetXaxis().SetTitle(
			"Number of sensor hits in box2")

			T0pesBOX2_histo_desc = "TimeFirstPESubtractedBOX2"
			self.T0pesBOX2_histo_name = self.alabel(T0pesBOX2_histo_desc)
			self.hman.h1(self.T0pesBOX2_histo_name, T0pesBOX2_histo_desc, 
			100, 0, 1000)
			self.hman.fetch(
			self.T0pesBOX2_histo_name).GetXaxis().SetTitle(
			"Time FirstPE Subtracted Box2 (ps)")

			TpesBOX2_histo_desc = "TimePESubtractedBox2"
			self.TpesBOX2_histo_name = self.alabel(TpesBOX2_histo_desc)
			self.hman.h1(self.TpesBOX2_histo_name, TpesBOX2_histo_desc, 
			100, 0, 1000)
			self.hman.fetch(
			self.TpesBOX2_histo_name).GetXaxis().SetTitle(
			"Time Inclusive PE Subtracted Box2 (ps)")

			DTpesBOX2_histo_desc = "TimeFirstPeMinusTimePeBox2"
			self.DTpesBOX2_histo_name = self.alabel(DTpesBOX2_histo_desc)
			self.hman.h1(self.DTpesBOX2_histo_name, DTpesBOX2_histo_desc, 
			100, 0, 1000)
			self.hman.fetch(
			self.DTpesBOX2_histo_name).GetXaxis().SetTitle(
			"Time Inclusive PE minus first PE Box2 (ps)")

			if self.QE <1:
				T0pesQeBOX2_histo_desc = "TimeFirstPESubtractedAfterQEBox2"
				self.T0pesQeBOX2_histo_name = self.alabel(T0pesQeBOX2_histo_desc)
				self.hman.h1(self.T0pesQeBOX2_histo_name, T0pesQeBOX2_histo_desc, 
				100, 0, 1000)
				self.hman.fetch(
				self.T0pesQeBOX2_histo_name).GetXaxis().SetTitle(
				"Time FirstPE Subtracted Box2 after QE (ps)")

				TpesQeBOX2_histo_desc = "TimePESubtractedAfterQEBox2"
				self.TpesQeBOX2_histo_name = self.alabel(TpesQeBOX2_histo_desc)
				self.hman.h1(self.TpesQeBOX2_histo_name, TpesQeBOX2_histo_desc, 
				100, 0, 1000)
				self.hman.fetch(
				self.TpesQeBOX2_histo_name).GetXaxis().SetTitle(
				"Time Inclusive PE Subtracted after QE Box2 (ps)")

			if self.SPTR > 0 or self.ASIC > 0 :
				T0pesQeSmearBOX2_histo_desc = "TimeFirstPESubtractedAfterQESmearBox2"
				self.T0pesQeSmearBOX2_histo_name = self.alabel(T0pesQeSmearBOX2_histo_desc)
				self.hman.h1(self.T0pesQeSmearBOX2_histo_name, T0pesQeSmearBOX2_histo_desc, 
				100, 0, 1000)
				self.hman.fetch(
				self.T0pesQeSmearBOX2_histo_name).GetXaxis().SetTitle(
				"Time FirstPE Subtracted Box2 after QE and smear (ps)")

				TpesQeSmearBOX2_histo_desc = "TimePESubtractedAfterQESmearBox2"
				self.TpesQeSmearBOX2_histo_name = self.alabel(TpesQeSmearBOX2_histo_desc)
				self.hman.h1(self.TpesQeSmearBOX2_histo_name, TpesQeSmearBOX2_histo_desc, 
				100, 0, 1000)
				self.hman.fetch(
				self.TpesQeSmearBOX2_histo_name).GetXaxis().SetTitle(
				"Time Inclusive PE Subtracted after QE and smear Box2 (ps)")

			
		
		
