"""
LXe
Properties of LXe
"""

from Xenon import *

class LXe:
  def __init__(self, wi=15.6*eV,ws=16.6*eV,lambdaScint=172*nm,rayleigh=36*mm,
      tau1=2.2*ns,tau2=27*ns,tau3=45*ns,
      rtau1=0.1,rtau2=0.2,rtau3=0.7,nUV=1.55,nBlue=1.4):
    
    self.Z = 54
    self.A = 131.29*g/mol
    self.T = 160*kelvin
    self.x0 = 8.48 * g/cm2
    self.rho = 3*g/cm3
    self.dedx=0.35*keV/micron  
    self.wi=wi
    self.ws=ws
    self.lambdaScint = lambdaScint
    self.rayleigh = rayleigh
    self.tau1=tau1
    self.tau2=tau2
    self.tau3=tau3
    self.rtau1=rtau1
    self.rtau2=rtau2
    self.rtau3=rtau3
    self.nUV = nUV
    self.nBlue=nBlue

  def AtomicNumber(self):
    """
    Xenon atomic number
    """
    return self.Z

  def AtomicMass(self):
    """
    Xenon atomic mass
    """
    return self.A

  def TemperatureAt1Bar(self):
    """
    LXe Temperature
    """
    return self.T
  def X0(self):
    """
    Xenon radiation length
    """
    return self.x0
  
  def RefractionIndexUV(self):
    return self.nUV

  def RefractionIndexBlue(self):
    return self.nBlue

  def Lifetime(self,i):
    """
    i ->(1,3) for the three lifetimes.
    """
    if i == 1: 
      return self.tau1
    elif i == 2: 
      return self.tau2
    elif i == 3: 
      return self.tau3
    else:
      print "index must be 1,2 or 3"
      sys.exit(0)

  def LifetimeRatio(self,i):
    """
    i ->(1,3) for the three lifetimes.
    """
    if i == 1: 
      return self.rtau1
    elif i == 2: 
      return self.rtau2
    elif i == 3: 
      return self.rtau3
    else:
      print "index must be 1,2 or 3"
      sys.exit(0)

  def Density(self):
    """
    Density of LXe 
    """
    return self.rho

  def Wi(self):
    """
    Energy needed to produce an ionization pair
    """
    return self.wi

  def Ws(self):
    """
    Energy needed to produce scintillation photons
    """
    return self.ws

  def ScintillationWavelength(self):
    """
    Scintillation wavelength 
    """
    return self.lambdaScint

  def Rayleigh(self):
    """
    Attenuation due to Rayleigh Scattering
    """
    return self.rayleigh

  def dEdX(self):
    return self.dedx
  
  def ComptonCrossSection(self,E):
    """
    Compton = Incoherent Scattering
    """
    return ScatterIncoherent(E)

  def PhotoelectricCrossSection(self,E):
    """
    Photoelectric Xs
    """
    return Photoelectric(E)

  def TotalCrossSection(self,E):
    """
    Total Xs
    """
    return TotalInteraction(E)

  def Attenuation(self,E,Z):
    """
    Attenuation of a beam of energy E in a thickness Z
    """
    return TransmittedBeam(E,Z,self.Density())
  
  def Efficiency(self,E,Z):
    """
    Fraction of gammas of energy E interacting in thickness E
    """
    return InteractionFraction(E,Z,self.Density())

  def GammaPathLength(self,E):
    """
    gamma path length in xenon 
    """
    xz = self.TotalCrossSection(E)*self.Density()
    return 1./xz

  def ScintillationPhotons(self,E):
    """
     Number of scintillation photons produced by a photon of energy E
    """
    return E/self.Ws()

  def SPhotonsAt511KeV(self,i):
    if i == 1: 
      return self.ScintillationPhotons(511*keV)*self.LifetimeRatio(1)
    elif i == 2: 
      return self.ScintillationPhotons(511*keV)*self.LifetimeRatio(2)
    elif i == 3: 
      return self.ScintillationPhotons(511*keV)*self.LifetimeRatio(3)
    else:
      print "index must be 1,2 or 3"
      sys.exit(0)

  def IonizationElectrons(self,E):
    """
     Number of ionization electrons produced by a photon of energy E
    """
    return E/self.Wi()

  def CostPerGram(self):
    """
    Cost per gram
    """
    return 1.0  #in euro


    
  def __str__(self):
    s= """
        Name = LXe  Z = %d  A = %7.4g g/mole 
        Temperature at atm pressure (1 bar) = %7.2f kelvin
        Density = %7.4g g/cm3 X0= %7.4g g/cm2 X1= %7.4g cm
        de/dx = %7.4g keV/cm Ws = %7.4g eV Wi = %7.4g eV
        Rayleigh Scattering = %7.2g cm
        Scintillation wavelength = %7.2f nm
        Refraction Index (UV) = %7.2f 
        Refraction Index (Blue) = %7.2f 
        Lifetimes:
        tau1 = %7.2f ns, ratio tau 1 = %7.2f
        tau2 = %7.2f ns, ratio tau 2 = %7.2f
        tau3 = %7.2f ns, ratio tau 3 = %7.2f 
        """%(self.AtomicNumber(),self.AtomicMass()/(g/mol),
          self.TemperatureAt1Bar()/kelvin,
          self.Density()/(g/cm3),
          self.X0()/(g/cm2),(self.X0()/self.Density())/cm,
          self.dEdX()/(keV/cm),self.Ws()/eV, self.Wi()/eV,
          self.Rayleigh()/cm, self.ScintillationWavelength()/nm, self.nUV, self.nBlue,
          self.tau1/ns, self.rtau1,self.tau2/ns,self.rtau2,
          self.tau3/ns, self.rtau3)

    return s 


if __name__ == '__main__':
    
    lxe = LXe()

    print lxe  

    print "Efficiency for 511 keV photons" 

    for z in drange(1., 11., 1.):
      print """
       z = %7.2g cm LXe eff = %7.2g 
      """%(z,
           lxe.Efficiency(511*keV,z*cm))

    print """Photoelectric fraction 
              at 511 keV photons = %7.2g"""%(
           lxe.PhotoelectricCrossSection(511*keV)/
           lxe.TotalCrossSection(511*keV))
    
    print """
       Gamma path lenght in LXe for 511 keV photons = %7.2g cm 
      """%(
           lxe.GammaPathLength(511*keV)/cm)

    
    print """
        Number of scintillation photons Ns (511 keV, LXe)= %7.2g 
        with tau1 = %7.2f ns lifetime: = %7.2f
        with tau2 = %7.2f ns lifetime: = %7.2f
        with tau3 = %7.2f ns lifetime: = %7.2f
      """%(
           lxe.ScintillationPhotons(511*keV),lxe.Lifetime(1),
           lxe.SPhotonsAt511KeV(1),
           lxe.Lifetime(2),
           lxe.SPhotonsAt511KeV(2),
           lxe.Lifetime(3),
           lxe.SPhotonsAt511KeV(3)
           )
           

        