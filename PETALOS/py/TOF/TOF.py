from Centella.physical_constants import *
import random as rnd
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

###########################################################
class SiPMHit(object):
	def __init__(self,hid,x,y,z,A,time):
		"""
		Describes a hit in a given SiPM
		id: id of SiPM
		x,y,z: position of SiPM
		A: ampliute (in pes)
		time: time stamp of hit
		"""
		self.id = hid
		self.hit = Point4D(x,y,z,time)
		self.A = A

	def T(self):
		return self.hit.T()

	def XYZ(self):
		return self.hit.XYZ()

	def __str__(self):
		s="""hit\n
			id = %s x = %7.2f mm, y = %7.2f mm z = %7.2f mm 
			t = %7.2f ps A = %7.2f pes
			"""%(self.id, self.hit.x, self.hit.y, self.hit.z,
				self.hit.t, self.A)

		return s
	

###########################################################
class TimeMap(object):
	def __init__(self,TimeMapBox1,TimeMapBox2):
		"""
		TimeMapBox1 and TimeMapBox2 are lists which describethe time-ordered zero-supressed
		waveform in a given SiPM. Each element of the list is a list 
		[hitId,(x,y,z,A,time)], where (x,y,z) are the hit position (SiPM coordinates)
		A is the amplitude and time is the hit time.  
		"""
    	
		self.tmb1=TimeMapBox1
		self.tmb2=TimeMapBox2
		self.nhb1 = len(self.tmb1)
		self.nhb2 = len(self.tmb2)

	def NumberOfHits(self):

		return (self.nhb1,self.nhb2)

	def timeHit(self,boxNumber=1,index=0,jitter=0):
		"""
		Returns the hit with index in box1 or box2 (boxNumber).
		index runs from 0 to the length of the TimeMap lists. Since the maps
		are ordered, '0' is the earliest time.  
		"""
		
		tmb = self.tmb1
		if boxNumber == 2:
			tmb = self.tmb2

		hit = tmb[index]
		hid = hit[0]
		values = hit[1]
		x = values[0]
		y = values[1]
		z = values[2]
		A = values[3]
		time = values[4]

		if jitter > 0:
			time += rnd.uniform(-1.,1.)*jitter

		#print "boxNumber =%d"%(boxNumber)
		#print "x,y,z,A,time",x,y,z,A,time

		return SiPMHit(hid,x,y,z,A,time)

	def __str__(self):
		s=' Time Map Box1\n'
		for hit in self.tmb1:
			hid = hit[0]
			values = hit[1]
			time = values[4]/ps
			s+="(ID = %s time =%7.2f ps); "%(hid,time)
		
		s+='\n Time Map Box2\n'
		for hit in self.tmb2:
			hid = hit[0]
			values = hit[1]
			time = values[4]/ps
			s+="(ID = %s time =%7.2f ps); "%(hid,time)

		return s


	