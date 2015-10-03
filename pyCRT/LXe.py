"""
Xenon
Properties of Xenon
"""

from Util import *
from abc import ABCMeta, abstractmethod
import sys

#Photon   ,Coherent,Incoher.,Photoel.,Nuclear ,Electron,Tot. w/ ,
#Energy   ,Scatter.,Scatter.,Absorb. ,Pr. Prd.,Pr. Prd.,Coherent,

Xe={
1.000E-03:[8.458E+00,4.417E-03,9.403E+03,0.000E+00,0.000E+00,9.412E+03],
1.072E-03:[8.393E+00,4.969E-03,8.133E+03,0.000E+00,0.000E+00,8.141E+03],
1.149E-03:[8.325E+00,5.582E-03,7.032E+03,0.000E+00,0.000E+00,7.040E+03],
1.149E-03:[8.325E+00,5.582E-03,7.334E+03,0.000E+00,0.000E+00,7.343E+03],
1.500E-03:[8.000E+00,8.522E-03,4.077E+03,0.000E+00,0.000E+00,4.085E+03],
2.000E-03:[7.481E+00,1.285E-02,2.080E+03,0.000E+00,0.000E+00,2.088E+03],
3.000E-03:[6.481E+00,2.123E-02,7.715E+02,0.000E+00,0.000E+00,7.780E+02],
4.000E-03:[5.637E+00,2.874E-02,3.730E+02,0.000E+00,0.000E+00,3.787E+02],
4.782E-03:[5.087E+00,3.383E-02,2.357E+02,0.000E+00,0.000E+00,2.408E+02],
4.782E-03:[5.087E+00,3.383E-02,6.890E+02,0.000E+00,0.000E+00,6.941E+02],
5.000E-03:[4.945E+00,3.515E-02,6.344E+02,0.000E+00,0.000E+00,6.393E+02],
5.104E-03:[4.880E+00,3.577E-02,5.995E+02,0.000E+00,0.000E+00,6.044E+02],
5.104E-03:[4.880E+00,3.577E-02,8.133E+02,0.000E+00,0.000E+00,8.182E+02],
5.275E-03:[4.775E+00,3.677E-02,7.515E+02,0.000E+00,0.000E+00,7.563E+02],
5.453E-03:[4.669E+00,3.777E-02,6.945E+02,0.000E+00,0.000E+00,6.992E+02],
5.453E-03:[4.669E+00,3.777E-02,8.018E+02,0.000E+00,0.000E+00,8.065E+02],
6.000E-03:[4.368E+00,4.071E-02,6.330E+02,0.000E+00,0.000E+00,6.374E+02],
8.000E-03:[3.464E+00,5.027E-02,2.998E+02,0.000E+00,0.000E+00,3.033E+02],
1.000E-02:[2.815E+00,5.848E-02,1.662E+02,0.000E+00,0.000E+00,1.691E+02],
1.500E-02:[1.852E+00,7.440E-02,5.550E+01,0.000E+00,0.000E+00,5.743E+01],
2.000E-02:[1.322E+00,8.486E-02,2.510E+01,0.000E+00,0.000E+00,2.651E+01],
3.000E-02:[7.545E-01,9.729E-02,8.078E+00,0.000E+00,0.000E+00,8.929E+00],
3.456E-02:[6.119E-01,1.008E-01,5.417E+00,0.000E+00,0.000E+00,6.130E+00],
3.456E-02:[6.119E-01,1.008E-01,3.244E+01,0.000E+00,0.000E+00,3.316E+01],
4.000E-02:[4.913E-01,1.038E-01,2.211E+01,0.000E+00,0.000E+00,2.270E+01],
5.000E-02:[3.492E-01,1.073E-01,1.227E+01,0.000E+00,0.000E+00,1.273E+01],
6.000E-02:[2.610E-01,1.091E-01,7.454E+00,0.000E+00,0.000E+00,7.824E+00],
8.000E-02:[1.608E-01,1.096E-01,3.363E+00,0.000E+00,0.000E+00,3.633E+00],
1.000E-01:[1.089E-01,1.081E-01,1.793E+00,0.000E+00,0.000E+00,2.010E+00],
1.500E-01:[5.298E-02,1.019E-01,5.651E-01,0.000E+00,0.000E+00,7.200E-01],
2.000E-01:[3.136E-02,9.555E-02,2.490E-01,0.000E+00,0.000E+00,3.759E-01],
3.000E-01:[1.464E-02,8.495E-02,8.009E-02,0.000E+00,0.000E+00,1.797E-01],
4.000E-01:[8.435E-03,7.688E-02,3.699E-02,0.000E+00,0.000E+00,1.223E-01],
5.000E-01:[5.472E-03,7.064E-02,2.088E-02,0.000E+00,0.000E+00,9.699E-02],
6.000E-01:[3.837E-03,6.559E-02,1.338E-02,0.000E+00,0.000E+00,8.281E-02],
8.000E-01:[2.181E-03,5.784E-02,6.940E-03,0.000E+00,0.000E+00,6.696E-02],
1.000E+00:[1.404E-03,5.211E-02,4.335E-03,0.000E+00,0.000E+00,5.785E-02],
1.022E+00:[1.345E-03,5.156E-02,4.126E-03,0.000E+00,0.000E+00,5.703E-02],
1.250E+00:[9.022E-04,4.665E-02,2.782E-03,1.914E-04,0.000E+00,5.052E-02],
1.500E+00:[6.279E-04,4.244E-02,1.991E-03,8.853E-04,0.000E+00,4.594E-02],
2.000E+00:[3.540E-04,3.625E-02,1.211E-03,2.971E-03,0.000E+00,4.078E-02],
2.044E+00:[3.390E-04,3.580E-02,1.169E-03,3.173E-03,0.000E+00,4.048E-02],
3.000E+00:[1.576E-04,2.854E-02,6.422E-04,7.458E-03,9.981E-06,3.681E-02],
4.000E+00:[8.871E-05,2.381E-02,4.250E-04,1.141E-02,4.070E-05,3.577E-02],
5.000E+00:[5.679E-05,2.057E-02,3.142E-04,1.481E-02,8.100E-05,3.583E-02],
6.000E+00:[3.944E-05,1.818E-02,2.478E-04,1.774E-02,1.243E-04,3.634E-02],
7.000E+00:[2.898E-05,1.634E-02,2.040E-04,2.035E-02,1.674E-04,3.709E-02],
8.000E+00:[2.219E-05,1.488E-02,1.731E-04,2.269E-02,2.091E-04,3.797E-02],
9.000E+00:[1.754E-05,1.367E-02,1.501E-04,2.482E-02,2.489E-04,3.891E-02],
1.000E+01:[1.420E-05,1.267E-02,1.324E-04,2.677E-02,2.868E-04,3.987E-02],
}

RHO={1:0.005395,2:0.010848,3:0.016361,4:0.02193,5:0.02757,6:0.03327,
7:0.03946,8:0.04488,9:0.05079,10:0.05678,
20:0.121127,
30:0.19682,40:0.289980,50:0.41501,60:0.62570}
WS={20:76.,60:25.6}

XENON={"Z":54,"A":131.29*g/mol,"X0":8.48 * g/cm2}

def selectE(E):
  """
  Select the key closest to the energy E
  """
  dE=999999.
  EE = -999999.

  #print "selectE: E = ", E/MeV
  for key in Xe.keys():
    diff=abs(E - key)
    if diff < dE:
      dE=diff
      EE=key 

  #print "EE =",EE
  return EE

def ScatterCoherent(E):
  """
  cross section due to coherent scattering
  """
  return Xe[selectE(E/MeV)][0]*cm2/g

def ScatterIncoherent(E):
  """
  cross section  due to incoherent scattering
  """
  return Xe[selectE(E/MeV)][1]*cm2/g

def Photoelectric(E):
  """
  cross section  due to photoelectric
  """
  return Xe[selectE(E/MeV)][2]*cm2/g

def PairNuclearFied(E):
  """
  cross section  due to pairs in nuclear field
  """
  return Xe[selectE(E/MeV)][3]*cm2/g
  
def PairElectronField(E):
  """
  cross section  due to pairs in electron field
  """
  return Xe[selectE(E/MeV)][4]*cm2/g
  
def TotalInteraction(E):
  """
   Total cross section 
    """ 
  return Xe[selectE(E/MeV)][5]*cm2/g  

def TransmittedBeam(E,Z,rho):
  """
  Transmitted beam of energy z through a thickness Z
  rho is the density
  """
  mu_over_rho =TotalInteraction(E)

  # print "mu_over_rho = %7.4g (cm2/g)"%(mu_over_rho/(cm2/g))
  # print "rho = %7.4g (g/cm3)"%(self.rho/(g/cm3))
  # print "z = %7.4g (cm)"%(Z/cm)
  # print "mu_over_rho*rho*z = %7.4g"%(mu_over_rho*self.rho*Z)
  # print "exp(-mu_over_rho*rho*z) = %7.4g"%(exp(-mu_over_rho*self.rho*Z))
  return exp(-mu_over_rho*Z*rho)

def InteractionFraction(E,Z,rho):
  """
  Fraction of gammas of energy E interacting in thickness E
  """
  return 1. - TransmittedBeam(E, Z,rho)
    
    

class AXenon:
  """
  Abstract class defining the interface functions for Xenon
  """
  __metaclass__ = ABCMeta
  
  
  @abstractmethod
  def AtomicNumber(self):
    """
    Xenon atomic number
    """
    pass

  @abstractmethod
  def AtomicMass(self):
    """
    Xenon atomic mass
    """
    pass
  @abstractmethod
  def Density(self):
    """
    Xenon density
    """
    pass

  @abstractmethod
  def Pressure(self):
    """
    Xenon Pressure
    """
    pass

  @abstractmethod
  def Temperature(self):
    """
    Xenon Temperature
    """
    pass

  def CostPerGram(self):
    """
    Cost per gram
    """
    pass
    
  @abstractmethod
  def Wi(self):
    """
    Energy needed to produce an ionization pair
    """
    pass

  @abstractmethod
  def Ws(self):
    """
    Energy needed to produce scintillation photons
    """
    pass
    
  @abstractmethod
  def X0(self):
    """
    radiation length
    """
    pass

  @abstractmethod
  def dEdX(self):
    """
    dE/dx
    """
    pass

  @abstractmethod
  def ComptonCrossSection(self):
    """
    Compton = Incoherent Scattering
    """
    pass

  @abstractmethod
  def PhotoelectricCrossSection(self):
    """
    Photoelectric Xs
    """
    pass

  @abstractmethod
  def TotalCrossSection(self):
    """
    Total Xs
    """
    pass

  @abstractmethod
  def Attenuation(self,E,Z):
    """
    Attenuation of a beam of energy E in a thickness Z
    """
    pass
  
  @abstractmethod
  def Efficiency(self,E,Z):
    """
    Fraction of gammas of energy E interacting in thickness E
    """
    pass

  @abstractmethod
  def GammaPathLength(self,E):
    """
    gamma path length in xenon 
    """
    pass

  @abstractmethod
  def ScintillationPhotons(self,E):
    """
     Number of scintillation photons produced by a photon of energy E
    """
    pass

  @abstractmethod
  def IonizationElectrons(self,E):
    """
     Number of ionization electrons produced by a photon of energy E
    """
    pass

    
class Xenon(AXenon):

  def AtomicNumber(self):
    """
    Xenon atomic number
    """
    return XENON["Z"]

  def AtomicMass(self):
    """
    Xenon atomic mass
    """
    return XENON["A"]

  def Density(self):
    """
    Xenon density
    """
    pass

  def Pressure(self):
    """
    Xenon Pressure
    """
    pass

  def Temperature(self):
    """
    Xenon Temperature
    """
    pass
    
  def Wi(self):
    """
    Energy needed to produce an ionization pair
    """
    pass

  def Ws(self):
    """
    Energy needed to produce scintillation photons
    """
    pass
    

  def X0(self):
    """
    Xenon radiation length
    """
    return XENON["X0"]

  
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



class LXe(Xenon):
  def __init__(self, wi=15.6*eV,ws=16.6*eV,lambdaAtt=36*mm,
      tau1=2.2*ns,tau2=27*ns,tau3=45*ns,
      rtau1=0.1,rtau2=0.2,rtau3=0.7,nuV=1.55,nBlue=1.4):
    
    
    self.rho = 3*g/cm3
    self.dedx=0.35*keV/micron  
    self.wi=wi
    self.ws=ws
    self.lambdaAtt = lambdaAtt
    self.tau1=tau1
    self.tau2=tau2
    self.tau3=tau3
    self.rtau1=rtau1
    self.rtau2=rtau2
    self.rtau3=rtau3
    self.nUV = nUV
    self.nBlue=nBlue
  
  def RefractionIndexUV(self):
    return nUV

  def RefractionIndexBlue(self):
    return nBlue

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
  def Pressure(self):
    """
    Notice that self.P is a scale variable but
    Pressure() returns pressure in physical units
    """
    return 1

  def Temperature(self):
    """
    Notice that self.P is a scale variable but
    Pressure() returns pressure in physical units
    """
    return 160*kelvin

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

  def LambdaAtt(self):
    return self.lambdaAtt

  def dEdX(self):
    return self.dedx

    
  def __str__(self):
    s= """
        Name = LXe  Z = %d  A = %7.4g g/mole 
        Density = %7.4g g/cm3 X0= %7.4g g/cm2 X1= %7.4g cm
        de/dx = %7.4g keV/cm Ws = %7.4g eV Wi = %7.4g eV
        LambdaAtt = %7.2g cm
        Refraction Index (UV) = %7.2f 
        Refraction Index (Blue) = %7.2f 
        Lifetimes:
        tau1 = %7.2f ns, ratio tau 1 = %7.2f
        tau2 = %7.2f ns, ratio tau 2 = %7.2f
        tau3 = %7.2f ns, ratio tau 3 = %7.2f 
        """%(self.AtomicNumber(),self.AtomicMass()/(g/mol),
          self.Density()/(g/cm3),
          self.X0()/(g/cm2),(self.X0()/self.Density())/cm,
          self.dEdX()/(keV/cm),self.Ws()/eV, self.Wi()/eV,
          self.LambdaAtt()/cm, self.nUV, self.nBlue,
          self.tau1/ns, self.rtau1,self.tau2/ns,self.rtau2,
          self.tau3/ns, self.rtau3)

    return s 


if __name__ == '__main__':
    
    lxe = LXe()

    print lxe  

    print "Efficiency for 140 keV photons" 

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

    print """Photoelectric fraction 
              at 511 keV photons = %7.2g"""%(
           lxe.ComptonCrossSection(511*keV)/
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
           lxe.ScintillationPhotons(511*keV)*lxe.LifetimeRatio(1),
           lxe.Lifetime(2),
           lxe.ScintillationPhotons(511*keV)*lxe.LifetimeRatio(2),
           lxe.Lifetime(3),
           lxe.ScintillationPhotons(511*keV)*lxe.LifetimeRatio(3)
           )
           

        