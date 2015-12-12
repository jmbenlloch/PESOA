import random as rnd
from LXe import *
from Geometry import *
from Util import *
from Centella.histoManager import *
from Centella.treeManager import *
from Centella.messenger import *
from system_of_units import *
from NumpyAlgebra import *
import array 
from ROOT import TF1


class PhotonGenerator:
	"""
	Generates a Xenon scintillation photon
	"""
	def __init__(self,lxe,level=0):

		self.m = Messenger(level)
		self.lxe = lxe
		self.rtau1=self.lxe.rtau1
   		self.rtau2=self.lxe.rtau2
		self.rtau3=self.lxe.rtau3		
	
		self.FTAU1 = TF1("FTAU1","exp(-x/[0])",0,10*self.lxe.tau1/ns)
		self.FTAU1.SetParameter(0,self.lxe.tau1/ns)

		self.FTAU2 = TF1("FTAU2","exp(-x/[0])",0,10*self.lxe.tau2/ns)
		self.FTAU2.SetParameter(0,self.lxe.tau2/ns)

		self.FTAU3 = TF1("FTAU3","exp(-x/[0])",0,10*self.lxe.tau3/ns)
		self.FTAU3.SetParameter(0,self.lxe.tau3/ns)

		self.FR = TF1("FR","exp(-x/[0])",0,10*self.lxe.tau3/ns)
		self.FTAU3.SetParameter(0,self.lxe.tau3/ns)

		self.m.log(1,'lxe = %s  '%(lxe))
		self.ntau1 = 0.
		self.ntau2 = 0.
		self.ntau3 = 0.
		self.ntot=0.
			
		
	def GeneratePhotonUV(self):
		"""
		Generates a VUV photon 
		"""
		t = 0
		tau=0
		self.ntot+=1
		test  = rnd.uniform(0.,1.)
		if test <= self.rtau1:
			t = self.FTAU1.GetRandom()*ns
			tau=1
			self.ntau1+=1
		elif test >self.rtau1 and test <= (self.rtau1+self.rtau2):
			t = self.FTAU2.GetRandom()*ns
			tau=2
			self.ntau2+=1
		else :
			t = self.FTAU3.GetRandom()*ns
			tau=3
			self.ntau3+=1

		self.m.log(2,'In GeneratePhotonUV: t0 = %7.2f ns tau = %d  '%(
			t/ns, tau))
		
		
		return (t,tau)

	def Statistics(self):
		return (self.ntau1/self.ntot,self.ntau2/self.ntot,self.ntau3/self.ntot)

	
		
	def __str__(self):
        
		s= """
        	lxe = %s

		"""%(self.lxe)
		return s



if __name__ == '__main__':
	nphotons = 200000
	nprint=100
	lxe = LXe()
	m = Messenger(1)

	hman =HistoManager() 
	pg = PhotonGenerator(lxe,level=1)

	hman.h1("tau1", "tau1", 
		100, 0, 10)
	hman.fetch("tau1").GetXaxis().SetTitle(
		"tau1 (ns)")

	hman.h1("tau2", "tau2", 
		100, 0, 150)
	hman.fetch("tau2").GetXaxis().SetTitle(
		"tau2 (ns)")

	hman.h1("tau3", "tau3", 
		100, 0, 250)
	hman.fetch("tau3").GetXaxis().SetTitle(
		"tau3 (ns)")

	hman.h1("tau23", "tau23", 
		100, 0, 250)
	hman.fetch("tau23").GetXaxis().SetTitle(
		"effective tau2-tau3 (ns)")

	hman.h1("tau12", "tau12", 
		100, 0, 150)
	hman.fetch("tau12").GetXaxis().SetTitle(
		"effective tau1-tau2 (ns)")

	hman.h1("tau123", "tau123", 
		100, 0, 250)
	hman.fetch("tau123").GetXaxis().SetTitle(
		"effective tau1-tau2-tau3 (ns)")

	for photons in xrange(nphotons):
		t,tau = pg.GeneratePhotonUV()
		if photons%nprint:
			print "number of photons generated = %d"%(photons)

		if tau == 1:
			hman.fill("tau1", t/ns)
			hman.fill("tau12", t/ns)
			hman.fill("tau123", t/ns)
		elif tau == 2:
			hman.fill("tau2", t/ns)
			hman.fill("tau12", t/ns)
			hman.fill("tau23", t/ns)
			hman.fill("tau123", t/ns)
		else:
			hman.fill("tau3", t/ns)
			hman.fill("tau123", t/ns)
			hman.fill("tau23", t/ns)

	print """
	fractions of events with lifetime tau1, tau2, tau3:
	ntau1 = %7.2f ntau2 = %7.2f ntau3 = %7.2f 
	"""%(pg.Statistics())

	pathFile = '/Users/jjgomezcadenas/Development/PETALO/WORK/histo/'
	fileName = 'LXeLifetime.root'
	hfile = pathFile+fileName
	hman.save(file_name=hfile)


	
	

