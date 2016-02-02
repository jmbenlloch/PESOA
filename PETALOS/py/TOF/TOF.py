from Util import *
import random as rnd
import sys
from Geometry import *

############################################################
def sortTimesSipm(hit):
	"""
	A helper function used to sort the hits. The hits are organises like this:
	(time, SipmHit): the function returns the time as key to sort the hit list
	"""
	return hit[0]

############################################################
def sortHits(hit):
	"""
	A helper function used to sort the hits. The hits are organises like this:
	(id, [x,y,z,A,t]): the function returns the time as key to sort the hit list
	"""
	hitVal = hit[1]
	return hitVal[4]

############################################################
def sortSiPmHits(hit):
	"""
	A helper function used to sort the SiPMhits. The hits are organises like this:
	(id, [x,y,z,A,t]): the function returns the time as key to sort the hit list
	"""
	sipm = hit[1]
	return sipm.W[0]

###########################################################
def SmearTime(time,SPTR,ctrASIC):
	"""
	Smears time to take into account the effect of SiPM and ASIC
	"""
	stime = time + rnd.gauss(0, SPTR)
	stime+= rnd.gauss(0, ctrASIC)
	return stime 


###########################################################
class SiPMHit(object):
	def __init__(self,hid,x,y,z,A,time,qE,sptr,ctrASIC,dT):
		"""
		Describes a hit in a given SiPM
		id: id of SiPM
		x,y,z: position of SiPM
		A: amplitude in pes for QE = 1. That is, number of VUV photons arriving to SiPM 
		time: time stamp of hit (first pe arriving to the SiPM, QE=1)
			  that is time stamp of the first VUV photon arriving to the SiPM. 
		W: a waveform containing the time stamp of all pes (after taking into account QE)
		   which are closer than DT to the first pe. 
		Q: the total charge in pes (after taking into account QE)
		qE is the quantum efficiency of the SiPM
		sptr is the single photon timr resolution of the SiPM, tipically 80 ps (rms)
		ctrASIC is the contribution ctr from the ASIC, tipically 30 ps (rms)
		dT is the time window (wrt first pes) to store pes in waveform 
		"""
		self.id = hid
		self.hit = Point4D(x,y,z,time)
		self.A = A
		self.W =[] #waveform
		self.Q = 0.
		self.dT = dT
		self.qE = qE
		self.sptr = sptr
		self.ctrASIC = ctrASIC

	def ID(self):
		return self.id

	def QE(self):
		return self.qE

	def DT(self):
		return self.dT

	def SPTR(self):
		return self.sptr

	def CtrASIC(self):
		return self.ctrASIC

	def TimeFirstUVPhoton(self):
		return self.hit.T()

	def NumberOfUVPhoton(self):
		return self.A

	def TimeFirstPE(self):
		if len(self.W) > 0:
			return self.W[0]
		else:
			return 0.

	def NumberOfPE(self):
		return self.Q

	def Waveform(self):
		return self.W

	def XYZ(self):
		return self.hit.XYZ()

	def __str__(self):
		
		s="""SiPMhit\n
			id = %s, QE = %7.2f DT = %7.2f ps 
			x = %7.2f mm, y = %7.2f mm z = %7.2f mm 
			number of VUV photons = %7.2f  time first UV photon = %7.2f ps
			number of PES = %7.2f  time first PE = %7.2f ps 
			"""%(self.ID(), self.QE(), self.DT()/ps, 
				self.hit.x/mm, self.hit.y/mm, self.hit.z/mm,
				self.NumberOfUVPhoton(), self.TimeFirstUVPhoton()/ps,
				self.NumberOfPE(), self.TimeFirstPE()/ps
				)
		s+="Waveform (ns)= %s"%(self.W)
		return s

		

###########################################################
class TimeMap(object):
	def __init__(self, numberOfBoxes=2, dtmax=1e+6):
		"""
		Class TimeMap organises the information in the boxes for TOF calculations

		vertexList: A list of interaction vertices for the gammas (one entry per box)
		each element of the vertexList is a Point4D which contains the true (x,y,z,t)
		of the interaction. 

		sipmMapList: a list of sipmMaps (one entry per box)
		each sipmMap is a time-ordered list [sipmID, sipmHit], where sipmID is the id
		of the SiPM and sipmHit is an instance of a class SiPMHit representing a hit
		in the SiPM.

		sipmMapDTList: same as sipmMapList, but only sipmHits with time within dtmax of 
		the first PE are included. 

		sipmMapTimeList: a list of sipmTimeMaps (one entry per box) 
		each sipmTime Map is a time-ordered list [time, sipmHit], where time is the time
		a PES and sipmHit is an instance of a class SiPMHit representing a hit
		in the SiPM. The list if built from the sipmMapDTList and contains an entry per each 
		PES whose time stamp is within dtmax of first pes.  


		"""
		self.numberOfBoxes = numberOfBoxes
		self.dtmax = dtmax
		self.vertexList =0
		self.siPmMapList =0
		self.siPmMapDTList =0
		self.siPmMapTimeList = 0
	
	def SetInteractionVertices(self,vertexList):
		"""
		Sets the vertex lists (The list contains one vertex per box)
		Each vertexList is a Point4d (x,y,z,t) describing the interaction point of the gamma 
		"""
		if (len(vertexList) != self.numberOfBoxes):
			print " error! length of vertex list =%d != number of boxes = %d"%(
				len(vertexList),self.numberOfBoxes)
			sys.exit()
		self.vertexList =vertexList

	def SetSiPmMaps(self,siPmMapList):
		"""
		Sets:
		self.siPmMapList --> one list per box. Each list is a time map. 
		a time map is a list [simpID, sipmObject], order by the time of this first PE
		of the sipmObject. 

		self.siPmMapDTList ---> one list per box, same structure, but the list only contains
		those sipm in an interval of dt wrt the first SiPM. 

		self.siPmMapTimeList  ---> one list per box. Each list is of the from
		[time, sipmObject] and has as many entries as pes in the box 
		"""
		if (len(siPmMapList) != self.numberOfBoxes):
			print " error! length of SiPM map list =%d != number of boxes = %d"%(
				len(siPmMapList),self.numberOfBoxes)
			sys.exit()

		self.siPmMapList = siPmMapList
		self.siPmMapDTList =self.setDTMaps(self.dtmax)
		self.siPmMapTimeList =self.setSiPmTimeMaps()

	
	def DTMAX(self):
		"""
		Return the DTMAX wrt first pes
		"""
		return self.dtmax

	def NumberOfBoxes(self):
		"""
		Return the number of boxes in the setup
		"""
		return self.numberOfBoxes

	def InteractionVertex(self, box=1):
		"""
		returns the interaction vertex in the box
		"""
		return self.vertexList[box-1]

	def SiPmMap(self, box=1):
		"""
		returns the SiPM map in the box
 		"""
		return self.siPmMapList[box-1]

	def SiPmMapDT(self, box=1):
		"""
		returns the SiPMDT map in the box
 		"""
		return self.siPmMapDTList[box-1]

	def SiPmTimeMap(self,box=1):
		"""
		returns the SiPM time map in box   
		"""
		return self.siPmMapTimeList[box-1]


	def SiPmHitId(self,box=1,index=0):
		"""
		Returns the SiPM hit id with index in box
		index runs from 0 to the number of SiPM in the box. Since the maps
		are ordered, '0' is the earliest time.  
		"""
		return self.SiPmMap(box)[index][0]

	def SiPmHit(self,box=1,index=0):
		"""
		Returns the SiPM hit with index in box
		index runs from 0 to the number of SiPM in the box. Since the maps
		are ordered, '0' is the earliest time.  
		"""
		return self.SiPmMap(box)[index][1]
	
	def NumberOfSiPmHits(self,box=1):
		return len(self.siPmMapList[box-1])

	def NumberOfSiPmHitsDT(self,box=1):
		return len(self.siPmMapDTList[box-1])

	def setDTMaps(self,dt):
		"""
		Creates a DT map list, considering only those SiPMs with a time stamp
		within dt of first SiPM hit.
		"""

		siPmMapDTList =[]
		for ibox in xrange(1,3):
			tSipmHit = self.SiPmHit(box=ibox,index=0).TimeFirstPE() 
			sipmMap = self.SiPmMap(box=ibox)
		
			hitMap =[]
			for elem in sipmMap:
				sipmid = elem[0]
				sipmhit = elem[1]
				if abs (sipmhit.TimeFirstPE() - tSipmHit) < dt:
					hitMaplist=[sipmid,sipmhit]
					hitMap.append(hitMaplist) 

			siPmMapDTList.append(hitMap)

		return siPmMapDTList


	def setSiPmTimeMaps(self):
		"""
		self.siPmMapTimeList  ---> one list per box. Each list is of the from
		[time, sipmObject] and has as many entries as pes in the box  
		"""
		siPmMapTimeList =[]

		for ibox in xrange(1,3): 
			sipmMap = self.SiPmMapDT(box=ibox)
			hitMap =[]
			for elem in sipmMap:
				sipmid = elem[0]
				sipmhit = elem[1]
				for tpes in sipmhit.W:
					hitMaplist=[tpes,sipmhit]
					hitMap.append(hitMaplist)
				
			siPmMapTimeList.append(hitMap)
		return siPmMapTimeList

	def displaySipPmMap(self,siPmMap):
		"""
		Display a siPmMap: ---> a time-ordered list [sipmID, sipmHit] 
		"""
		
		s='['
		for elem in siPmMap:
			hid = elem[0]
			sipm = elem[1]
			s+='(id=%d,time=%7.2f) '%(hid,sipm.TimeFirstPE())
		s+=']'
		
		return s

	def displaySiPmTimeMap(self,timeMap):
		"""
		Display a SipmTimeMap --> a time-ordered list [time, sipmHit]
		"""
		
		s='['
		for elem in timeMap:
			time = elem[0]
			sipm = elem[1]
			s+='(time=%7.2f,id=%d ) '%(time,sipm.ID())
		s+=']'
		"""
		returns the SiPM time map in box   
		"""
		return s

	def __str__(self):
		
		s=' Time Map: Number of boxes =%d, dtmax =%7.2f '%(
			self.NumberOfBoxes(), self.dtmax)
		for nb in xrange(0,self.NumberOfBoxes()):
			s+='\n-->box %d'%(nb+1)
			s+='\nlength of SiPmMap = %d, '%(len(self.SiPmMap(nb)))
			s+='length of SiPmMapDT = %d, '%(len(self.SiPmMapDT(nb)))
			s+='length of SiPmTimeMap = %d'%(len(self.SiPmTimeMap(nb)))
			s+='\n+++SiPmMap (DT)  \n '
			s+=self.displaySipPmMap(self.SiPmMapDT(box=nb))
			s+='\n ---SiPmTimeMap (DT)\n '
			s+=self.displaySiPmTimeMap(self.SiPmTimeMap(box=nb))
		
		return s
	
