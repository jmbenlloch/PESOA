"""
LXe
Properties of LXe
"""

from Xenon import *

#energy in eV
LXeRefractionIndex =[
[6.4, 1.58587, 0.0964027],
[6.6, 1.61513, 0.508607],
[6.8, 1.6505, 1.33957],
[7, 1.69447, 1.69005],
[7.2, 1.75124, 1.02138],
[7.4, 1.82865, 0.295683],
[7.6, 1.94333, 0.098714]
]

############################################################
def sortRI(elem):
  """
  A helper function used to sort the hits. The hits are organises like this:
  (id, [x,y,z,A,t]): the function returns the time as key to sort the hit list
  """
  return elem[0]

class LXe:
  def __init__(self, wi=15.6*eV,ws=16.6*eV,lambdaScint=172*nm,rayleigh=36*mm,
      tau1=2.2*ns,tau2=27*ns,tau3=45*ns,
      rtau1=0.065,rtau2=0.935,rtau3=0.0,nUV=1.70,nBlue=1.4):
    
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
    lxri =[]  #transform to nm
    for elem in LXeRefractionIndex:
      ene = elem[0]
      n = elem[1]
      f = elem[2]
      x =[(1240./ene)*nm,n,f]
      #print x[0]/nm
      lxri.append(x)

    self.LXRI = sorted(lxri, key=sortRI)
    #print self.LXRI

    l,n = self.AverageLamdaAndRI()
    self.AverageLamda = l
    self.AverageRefractionIndexUV = n

  def AverageLamdaAndRI(self):
    """
    Returns the average lamda and refraction index
    """
    l=0.
    n=0.
    w=0.
    for elem in self.LXRI:
      l+=elem[0]*elem[2]
      n+=elem[1]*elem[2]
      w+=elem[2]
    return (l/w,n/w)



  def RefractionIndex(self,lamda):
    """
    returns the refraction index
    """

    if lamda < self.LXRI[0][0]:
      return self.LXRI[0][1] 
    elif lamda > self.LXRI[6][0]:
      return self.LXRI[6][1]
    else:
      for i in xrange(len(self.LXRI)-1):
        elem = self.LXRI[i]
        x0 = elem[0]
        y0 = elem[1]
        elem = self.LXRI[i+1]
        x1 = elem[0]
        y1 = elem[1]
        if lamda >= x0 and lamda < x1:
          break
      return lin(lamda,x0,y0,x1,y1)

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
    return self.AverageRefractionIndexUV

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

    for l in drange(150*nm,220*nm,5*nm):
      print """
      for lamda = %7.2f nm (%7.2f eV) n = %7.2f
      """%(l/nm, 1240./(l/nm), lxe.RefractionIndex(l))

    l,n = lxe.AverageLamdaAndRI()
    print """
    Average lamda = %7.2f nm ; average n = %7.2f
    """%(l/nm, n)
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
           

        