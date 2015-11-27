from Centella.AAlgo import AAlgo
from Centella.physical_constants import *
from Particles import *

"""
This algorithm classifies the interactions in the box(es) 
according to the following types:
i) Photoelectric
ii) Compton
Parameters:
   pdfMinE -> PDF Minimum  Energy
   pdfMaxE -> PDF Maximum  Energy
   roiMinE -> ROI Minimum  Energy
   roiMaxE -> ROI Maximum  Energy
"""

class CEvent(AAlgo):

	############################################################
	def __init__(self, param=False, level=1, label="", **kargs):

		"""
		CEvent Algorithm
		"""
		#self.m.log(1, 'Constructor()')

		### GENERAL STUFF
		self.name = 'CEvent'
		
		AAlgo.__init__(self, param, level, self.name, 0, label, kargs)

    ### PARAMETERS
    # Box coordinates

  		self.boxCoord =[]
  		print self.vdoubles
		
		for i in range(1,9):
			print i
			blb="BoxV"
			blb +=str(i)
			print blb
			self.m.log(1, "loading parameter", blb)
			coord  = self.vdoubles[blb]
			try:
				self.boxCoord.append(coord);
			except KeyError:
				self.m.log(1, "WARNING!! Parameter %s not defined."%(blb))
				exit(0)
  
		self.m.log(1, "Box coordinates", self.boxCoord)
		

   

	############################################################		
	def initialize(self):

		self.m.log(1, 'Initialize()')
		
		### Defining histos
		# Event energy histogram

		self.NumberOfParticles_histo_desc = "NumberOfParticles"
		self.NumberOfCompton_histo_desc = "NumberOfParticlesInCompton"
		self.NumberOfPhoto_histo_desc = "NumberOfParticlesInPhoto"
		self.NumberOfParticles_histo_name = self.alabel(self.NumberOfParticles_histo_desc)
		self.NumberOfCompton_histo_name = self.alabel(self.NumberOfCompton_histo_desc)
		self.NumberOfPhoto_histo_name = self.alabel(self.NumberOfPhoto_histo_desc)

		self.Edep_histo_desc = "Edep"
		self.EdepCompton_histo_desc = "EdepCompton"
		self.EdepPhoto_histo_desc = "EdepPhoto"
		
		self.Edep_histo_name = self.alabel(self.Edep_histo_desc)
		self.EdepCompton_histo_name = self.alabel(self.EdepCompton_histo_desc)
		self.EdepPhoto_histo_name = self.alabel(self.EdepPhoto_histo_desc)
		
		self.m.log(2, "Booking histograms ")
		
		self.hman.h1(self.NumberOfParticles_histo_name, self.NumberOfParticles_histo_desc, 
			20,0.,20.)
		self.hman.fetch(self.NumberOfParticles_histo_name).GetXaxis().SetTitle("Number of particles")

		self.hman.h1(self.NumberOfCompton_histo_name, self.NumberOfCompton_histo_desc, 
			20,0.,20.)
		self.hman.fetch(self.NumberOfCompton_histo_name).GetXaxis().SetTitle("Particles in Compton")

		self.hman.h1(self.NumberOfPhoto_histo_name, self.NumberOfPhoto_histo_desc, 
			20,0.,20.)
		self.hman.fetch(self.NumberOfPhoto_histo_name).GetXaxis().SetTitle("Particles in Photo")

		self.hman.h1(self.Edep_histo_name, self.Edep_histo_desc, 
			100,0.,600.)
		self.hman.fetch(self.Edep_histo_name).GetXaxis().SetTitle("Energy deposited (keV)")

		self.hman.h1(self.EdepCompton_histo_name, self.EdepCompton_histo_desc, 
			20,0.,600.)
		self.hman.fetch(self.EdepCompton_histo_name).GetXaxis().SetTitle("Energy deposited Compton (keV)")

		self.hman.h1(self.EdepPhoto_histo_name, self.EdepPhoto_histo_desc, 
			20,0.,600.)
		self.hman.fetch(self.EdepPhoto_histo_name).GetXaxis().SetTitle("Energy deposited Compton (keV)")

    ### Counters:
		self.numInputEvents = 0
		self.numOutputEvents = 0
		self.compton = 0
		self.photo = 0

		return



	############################################################
	def execute(self, event=""):

		self.m.log(2, 'Execute()')		
	
		self.numInputEvents += 1   

		evtEnergy = event.GetEnergy()
		numberOfParticles = event.GetMCParticles().size()

		self.m.log(3, ' number of particles = %d '%(numberOfParticles))
		self.hman.fill(self.NumberOfParticles_histo_name, numberOfParticles)

		primaryParticles = PrimaryParticles(event)
		self.m.log(3, ' number of primary Particles =%d '%(len(primaryParticles)))
		
		for pparticle in primaryParticles:
			self.m.log(3, ' +++primary particle+++')
			self.particleInfo(pparticle)

			secondaryParticles = SecondaryParticles(pparticle)
			interactionType = ClassifyInteraction(secondaryParticles[0])
			self.m.log(3, ' interaction type = %s '%(interactionType))

			Edep =0
			for sparticle in secondaryParticles:
				self.m.log(3, '--secondary particle--')
				self.particleInfo(sparticle)
				ti,tf = particleKineticEnergy(sparticle)
				Edep += ti 

			self.m.log(3, ' Edep =%7.2f keV '%(Edep/keV))

			self.hman.fill(self.NumberOfParticles_histo_name, len(secondaryParticles))
			self.hman.fill(self.Edep_histo_name, Edep/keV)
			if interactionType == "phot":
				self.hman.fill(self.NumberOfPhoto_histo_name, len(secondaryParticles))
				self.hman.fill(self.EdepPhoto_histo_name, Edep/keV)
			else:
				self.hman.fill(self.NumberOfCompton_histo_name, len(secondaryParticles))
				self.hman.fill(self.EdepCompton_histo_name, Edep/keV)
		if interactionType == "phot":
			self.photo +=1
		else:
			self.compton+=1
			return False

		self.numOutputEvents += 1 
		return True

		

	############################################################
	def finalize(self):

		self.m.log(1, 'Finalize()')

		self.m.log(1, 'Input  Events: ', self.numInputEvents)
		self.m.log(1, 'Output Events: ', self.numOutputEvents)
		self.m.log(1, 'Photoelectric Events: ', self.photo)
		self.m.log(1, 'Compton Events: ', self.compton)

		self.logman["USER"].ints[self.alabel("InputEvents")] = self.numInputEvents
		self.logman["USER"].ints[self.alabel("OutputEvents")] = self.numOutputEvents
		self.logman["USER"].ints[self.alabel("ComptonEvents")] = self.compton
		self.logman["USER"].ints[self.alabel("PhotoelectricEvents")] = self.photo

		return

	############################################################
	def particleInfo(self, particle):

		self.m.log(3,'name = %s t =%7.2f ns '%(
			particleName(particle), particleTime(particle)/ns))
		
		ei,ef = particleEnergy(particle)
		ti,tf = particleKineticEnergy(particle)
		self.m.log(3,'Ei = %7.2f keV ef = %7.2f keV Ti = %7.2f keV Tf= %7.2f keV '%(
			ei/keV,ef/keV,ti/keV,tf/keV)) 

		x,y,z = particleInitialVtx(particle)
		self.m.log(3, ' Initial vertex: x =%7.2f mm, y =%7.2f mm, z =%7.2f mm '%(
			x/mm,y/mm,z/mm))
		x,y,z = particleFinalVtx(particle)
		self.m.log(3, ' Final vertex: x =%7.2f mm, y =%7.2f mm, z =%7.2f mm '%(
			x/mm,y/mm,z/mm))

	