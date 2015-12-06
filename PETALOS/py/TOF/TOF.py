from Centella.physical_constants import *
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

	def timeHit(self,boxNumber=1,index=0):
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

		#print "boxNumber =%d"%(boxNumber)
		#print "x,y,z,A,time",x,y,z,A,time

		return (x,y,z,A,time)

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


	