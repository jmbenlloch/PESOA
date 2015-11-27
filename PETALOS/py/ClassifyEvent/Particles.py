from math import sqrt
############################################################
def PrimaryParticles(event):
	"""
	Finds the primary particles in the event
	(and returns them in a list sorted by time)
	"""
	
	primaryParticles =[]
	for particle in event.GetMCParticles():
		if particle.IsPrimary():
			
			#print '+++++found primary particle, name = ',particleName(particle)
			#particleInfo(particle)
			primaryParticles.append(particle)

	return sorted(primaryParticles, key=particleTime)

############################################################
def SecondaryParticles(primaryParticle):
	"""
	Finds the particles associated to a primaryParticle 
	(and returns them in a list sorted by time)
	"""
	
	Particles =[]

	for particle in primaryParticle.GetDaughters():
		#print '-----found secondary partic, name = ',particleName(particle) 
		#particleInfo(particle)
		Particles.append(particle)
	return sorted(Particles, key=particleTime)


############################################################
def ClassifyInteraction(particle):
	#print 'In ClassifyInteraction: process ', particle.GetCreatorProc()
	return particle.GetCreatorProc() 

############################################################
def particleTime(particle):
	return particle.GetInitialVtx4D().GetT()

############################################################
def particleEnergy(particle):
	return (particle.GetInitialMom().GetE(), particle.GetFinalMom().GetE())

############################################################
def particleKineticEnergy(particle):
	# print "particle energy =",particleEnergy(particle)
	# print "particle mass =",particleMass(particle)
	# print "kinetic energy =",particleEnergy(particle) - particleMass(particle)
	ei, ef = particleEnergy(particle)
	return (ei - particleMass(particle), ef - particleMass(particle))

############################################################
def particleName(particle):
	pdg = particle.GetPDG()
	if  pdg == 11:
  		return "e-"
  	elif pdg == -11 :
  		return "e+"
  	elif pdg == 22 :
  		return "gamma"
  	else :
  		return "unknown"
############################################################
def particleMass(particle):
	pdg = particle.GetPDG()
	if  pdg == 11:
  		return 0.510998902
  	elif pdg == -11 :
  		return 0.510998902
  	elif pdg == 22 :
  		return 0.
  	else :
  		return 0.


############################################################
def particleInitialVtx(particle):
	return (particle.GetInitialVtx().x(), particle.GetInitialVtx().y(),
	particle.GetInitialVtx().z())

############################################################
def particleFinalVtx(particle):
	return (particle.GetFinalVtx().x(), particle.GetFinalVtx().y(),
	particle.GetFinalVtx().z())
	
############################################################
def particleInfo(particle):
	print "***Particle Info***"
	ei,ef = particleEnergy(particle)
	ti,tf = particleKineticEnergy(particle)
	print " name = %s, mass = %7.2f, Ei = %7.2f Ef = %7.2f Ti = %7.2f Tf = %7.2f"%(
		particleName(particle), particleMass(particle),ei,ef,ti,tf)
	print 'time = %7.2f '%(particleTime(particle))		
	x,y,z = particleInitialVtx(particle)
	print 'initial vertex x = %7.2f y = %7.2f z = %7.2f '%(x,y,z)
	x,y,z = particleFinalVtx(particle)		
	print 'final vertex x = %7.2f y = %7.2f z = %7.2f '%(x,y,z)

	
