from Centella.physical_constants import *
import random as rnd
import sys
from Geometry import *
ps = picosecond

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
class SiPMHit(object):
	def __init__(self,hid,x,y,z,A,time,qE,dT):
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
		dT is the time window (wrt first pes) to store pes in waveform 
		"""
		self.id = hid
		self.hit = Point4D(x,y,z,time)
		self.A = A
		self.W =[] #waveform
		self.Q = 0.
		self.dT = dT
		self.qE = qE

	def ID(self):
		return self.id

	def QE(self):
		return self.qE

	def DT(self):
		return self.dT

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
	def __init__(self, numberOfBoxes=2):
		"""
		Class TimeMap organises the information in the boxes for TOF calculations

		vertexList: A list of interaction vertices for the gammas (one entry per box)
		each element of the vertexList is a Point4D which contains the true (x,y,z,t)
		of the interaction. 

		sipmMapList: a list of sipmMaps (one entry per box)
		each sipmMap is a time-ordered list [sipmID, sipmHit], where sipmID is the id
		of the SiPM and sipmHit is an instance of a class SiPMHit representing a hit
		in the SiPM. 
		"""
		self.numberOfBoxes = numberOfBoxes
		self.siPmMapList =0
		self.vertexList =0

	def SetSiPmMaps(self,siPmMapList):
		"""
		Sets the SiPM map lists (The list contains one SiPmMap per box)
		Each SiPmMap is a list [sipmId, sipmInstance]  
		"""
		if (len(siPmMapList) != self.numberOfBoxes):
			print " error! length of SiPM map list =%d != number of boxes = %d"%(
				len(siPmMapList),self.numberOfBoxes)
			sys.exit()

		self.siPmMapList = siPmMapList

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

	def SiPmHitId(self,box=1,index=0):
		"""
		Returns the SiPM hit with index in box
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
		
		# if jitter > 0:

		# 	time += rnd.gauss(0, jitter)

		#print "boxNumber =%d"%(boxNumber)
		#print "x,y,z,A,time",x,y,z,A,time


	def __str__(self):
		
		s=' Time Map: Number of boxes =%d '%(self.NumberOfBoxes())
		for nb in xrange(0,self.NumberOfBoxes()):
			s+=' box %d '%(nb)
			s+=' length of SiPmMap = %d'%(len(self.SiPmMap(nb)))
			for indx in xrange(0,len(self.SiPmMap(nb))):
				s+='[%s]'%(self.SiPmHit(box=nb,index=indx))
		return s
	