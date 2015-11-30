from Centella.AAlgo import AAlgo
from Centella.physical_constants import *
from Particles import *
from Util import *
from Geometry import *

"""
This algorithm classifies the interactions in the box(es) 
according to the following types:
i) Photoelectric
ii) Compton
Parameters:
   box coordinates
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

  		#print self.vdoubles

  		boxCoord1 =self.loadCoord("Box1V")
  		boxCoord2 =self.loadCoord("Box2V")

		self.box1 = Box(boxCoord1)
		self.box2 = Box(boxCoord2)
		
		self.m.log(2, "Box1 --", self.box1)
		self.m.log(2, "Box2 ---", self.box2)

		#wait()
   

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
		
		self.hman.h1(self.NumberOfParticles_histo_name, 
			self.NumberOfParticles_histo_desc, 
			20,0.,20.)
		self.hman.fetch(
			self.NumberOfParticles_histo_name).GetXaxis().SetTitle(
			"Number of particles")

		self.hman.h1(self.NumberOfCompton_histo_name, 
			self.NumberOfCompton_histo_desc, 
			20,0.,20.)
		self.hman.fetch(
			self.NumberOfCompton_histo_name).GetXaxis().SetTitle(
			"Particles in Compton")

		self.hman.h1(self.NumberOfPhoto_histo_name, 
			self.NumberOfPhoto_histo_desc, 
			20,0.,20.)
		self.hman.fetch(
			self.NumberOfPhoto_histo_name).GetXaxis().SetTitle(
			"Particles in Photo")

		self.hman.h1(self.Edep_histo_name, 
			self.Edep_histo_desc, 
			100,0.,600.)
		self.hman.fetch(
			self.Edep_histo_name).GetXaxis().SetTitle(
			"Energy deposited (keV)")

		self.hman.h1(self.EdepCompton_histo_name, 
			self.EdepCompton_histo_desc, 
			20,0.,600.)
		self.hman.fetch(
			self.EdepCompton_histo_name).GetXaxis().SetTitle(
			"Energy deposited Compton (keV)")

		self.hman.h1(self.EdepPhoto_histo_name, 
			self.EdepPhoto_histo_desc, 
			20,0.,600.)
		self.hman.fetch(
			self.EdepPhoto_histo_name).GetXaxis().SetTitle(
			"Energy deposited Compton (keV)")

		self.XYZBox1_histo_desc = "XYZBox1"
		self.XYZBox1_histo_name = self.alabel(self.XYZBox1_histo_desc)
		self.hman.h3(self.XYZBox1_histo_name, self.XYZBox1_histo_desc, 
			10, self.box1.xmin, self.box1.xmax,
           	10, self.box1.ymin, self.box1.ymax,
           	10, self.box1.zmin, self.box1.zmax)
		self.hman.fetch(
			self.XYZBox1_histo_name).GetXaxis().SetTitle(
			"XYZ interaction box 1 (mm)")


		self.XYZBox2_histo_desc = "XYZBox2"
		self.XYZBox2_histo_name = self.alabel(self.XYZBox2_histo_desc)
		self.hman.h3(self.XYZBox2_histo_name, self.XYZBox2_histo_desc, 
			10, self.box2.xmin, self.box2.xmax,
           	10, self.box2.ymin, self.box2.ymax,
           	10, self.box2.zmin, self.box2.zmax)
		self.hman.fetch(
			self.XYZBox2_histo_name).GetXaxis().SetTitle(
			"XYZ interaction box 2 (mm)")

    ### Counters:
		self.numInputEvents = 0
		self.numOutputEvents = 0
		self.compton = 0
		self.photo = 0
		self.ngbox0 = 0 # number of gammas not in box1 or box2
		self.ngbox1 = 0 # number of gammas in box1
		self.ngbox2 = 0 # number of gammas in box2

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
		
		npb0 = 0
		npe = 0
		nco = 0
		interactionType = "none"

		for pparticle in primaryParticles:
			self.m.log(3, ' +++primary particle+++')
			self.particleInfo(3,pparticle)
			ei,ef = particleEnergy(pparticle)

			if ei/keV != 511 : #loop over high energy gamma
				continue 

			box,x,y,z = self.primaryParticleFiducial(pparticle)
			if box == 1: 
				self.ngbox1+=1
			elif box == 2: 
				self.ngbox2+=1
			else:
				self.ngbox0+=1
				npb0+=1

			if box == 1:
				self.hman.fill(self.XYZBox1_histo_name, 
					x/mm,y/mm,z/mm)
			elif box == 2:
				self.hman.fill(self.XYZBox2_histo_name, 
					x/mm,y/mm,z/mm)
			else:
				continue

			secondaryParticles = SecondaryParticles(pparticle)
			
			if len(secondaryParticles) == 0: # loop over particles with no sec
				self.m.log(3, ' !!No secondaries!!')
				continue

			interactionType = ClassifyInteraction(secondaryParticles[0])
			self.m.log(3, ' interaction type = %s '%(interactionType))

			Edep =0
			for sparticle in secondaryParticles:
				self.m.log(4, '--secondary particle--')
				self.particleInfo(4,sparticle)
				ti,tf = particleKineticEnergy(sparticle)
				Edep += ti 

			self.m.log(3, ' Edep =%7.2f keV '%(Edep/keV))


			self.hman.fill(self.NumberOfParticles_histo_name, len(secondaryParticles))
			self.hman.fill(self.Edep_histo_name, Edep/keV)
			

			if interactionType == "phot":
				self.hman.fill(self.NumberOfPhoto_histo_name, len(secondaryParticles))
				self.hman.fill(self.EdepPhoto_histo_name, Edep/keV)
				npe+=1
				self.photo +=1
			elif interactionType == "compt":
				self.hman.fill(self.NumberOfCompton_histo_name, len(secondaryParticles))
				self.hman.fill(self.EdepCompton_histo_name, Edep/keV)
				self.compton+=1
				nco+=1
		
		self.m.log(2, ' number of gammas out of box1 and box2 =%d  '%(npb0))
		self.m.log(2, ' number of gammas that do photoelectric =%d  '%(npe))
		self.m.log(2, ' number of gammas that do compton =%d  '%(nco))
		
		#wait()
		if npb0 ==0 and npe == 2: # both gammas interact and are photo
			return True
			self.numOutputEvents += 1
		else:
			return False

		

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
	def particleInfo(self, lvl, particle):
		self.m.log(lvl,' ++ParticleInfo++')
		self.m.log(lvl,'name = %s t =%7.2f ns '%(
			particleName(particle), particleTime(particle)/ns))
		
		ei,ef = particleEnergy(particle)
		ti,tf = particleKineticEnergy(particle)
		self.m.log(lvl,'Ei = %7.2f keV ef = %7.2f keV Ti = %7.2f keV Tf= %7.2f keV '%(
			ei/keV,ef/keV,ti/keV,tf/keV)) 

		x,y,z = particleInitialVtx(particle)
		self.m.log(lvl, ' Initial vertex: x =%7.2f mm, y =%7.2f mm, z =%7.2f mm '%(
			x/mm,y/mm,z/mm))
		x,y,z = particleFinalVtx(particle)
		self.m.log(lvl, ' Final vertex: x =%7.2f mm, y =%7.2f mm, z =%7.2f mm '%(
			x/mm,y/mm,z/mm))

	############################################################
	def primaryParticleFiducial(self, particle):

		self.m.log(3,' ++primaryParticleFiducial++')
		x,y,z = particleFinalVtx(particle)
		self.m.log(3, ' Final vertex: x =%7.2f mm, y =%7.2f mm, z =%7.2f mm '%(
			x/mm,y/mm,z/mm))

		box = 0 # by default out of box1 and box2
		fid = self.box1.Active([x,y,z])
		if fid:
			self.m.log(3,'gamma found in box1')
			box = 1
		else:
			fid = self.box2.Active([x,y,z])
			if fid:
				self.m.log(3,'gamma found in box2')
				box = 2

		if box == 0:
			self.m.log(3,'gamma not found in box1 and in box2')

		return (box,x,y,z)

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

	