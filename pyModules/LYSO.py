"""
LYSO
Properties of LYSO
"""
from Centella.physical_constants import *
import sys
from Util import *


nm = nanometer 
mol = mole
micron = micrometer

LysoScintSpectrum =[
[2.2543, 0.0643],
[2.4797, 0.1929],
[2.6436, 0.4],
[2.8700, 1.],
[3.0996, 0.4071]
]
def sortRI(elem):
  """
  A helper function used to sort the hits. The hits are organises like this:
  (id, [x,y,z,A,t]): the function returns the time as key to sort the hit list
  """
  return elem[0]

class LYSO:
    def __init__(self,Z=54,rho=7.3*g/cm3, n=1.82, X0=1.16*cm, LambdaAtt = 0.87*(1./cm),
        LambdaPeak=420*nm, tau = 50*ns, PhotoFraction = 0.3, Nphot = 15000):
    
        """
        Represents lyso
        """  
        self.x0 = X0  
        self.Z = Z
        self.rho = rho
        self.tau=tau
        self.n = n
        self.mu=LambdaAtt
        self.lambdaScint = LambdaPeak
        self.tau = tau
        self.photoF = PhotoFraction
        self.Nphot = Nphot

        lysct =[]  #transform to nm
        for elem in LysoScintSpectrum:
            ene = elem[0]
            w =elem[1]
            x =[(1240./ene)*nm,w]
      #print x[0]/nm
            lysct.append(x)

        self.LYSC = sorted(lysct, key=sortRI)
    
        print "scintillation spectrum"
        for elem in self.LYSC:
            print " lambda = %7.2f nm w= %7.2g"%(elem[0]/nm,elem[1])

        l = self.AverageLamda()
        
    
    def AverageLamda(self):
        """
        Returns the average lamda 
        """
        l=0.
        w=0.
        for elem in self.LYSC:
            l+=elem[0]*elem[1]
            w+=elem[1]
        return (l/w)

    def EffectiveAtomicNumber(self):
        """
        EffectiveAtomicNumber
        """
        return self.Z

    def X0(self):
        """
        EffectiveAtomicNumber
        """
        return self.x0

    def RefractionIndex(self):
        """
        returns the refraction index
        """
        return self.n

    def Lifetime(self):
        """
        returns the lifetime
        """
        return self.tau 

    def Density(self):
        """
        Density 
        """
        return self.rho

    def ScintillationWavelength(self):
        """
        Scintillation wavelength 
        """
        return self.lambdaScint

    def PhotoelectricFraction(self):
        """
        Photoelectric Xs
        """
        return self.photoF

    def Attenuation(self,Z):
        """
        Attenuation of a beam of energy E 511 keV in a thickness Z
        """
        return exp(-Z*self.mu)
  
    def Efficiency(self,Z):
        """
        Fraction of gammas of energy E 511 keV interacting in thickness Z
        """
        return 1. - self.Attenuation(Z)

    def ScintillationPhotons(self):
        """
        Number of scintillation photons produced by a photon of energy 511 keV
        """
        return self.Nphot

    def __str__(self):
        s= """
            Name = LYSO  Z = %d  
            Density = %7.4g g/cm3 X0= %7.4g g/cm2 
            Scintillation wavelength = %7.2f nm
            Refraction Index (Blue) = %7.2f 
            Lifetime = %7.2f ns
            ScintillationPhotons =  %7.2f
            Attenuation in 1 cm =  %7.2f
            PhotoelectricFraction = %7.2f
            """%(self.EffectiveAtomicNumber(),
            self.Density()/(g/cm3),
            self.X0()/cm,
            self.ScintillationWavelength()/nm, self.RefractionIndex(),
            self.Lifetime(), self.ScintillationPhotons(),
            self.Attenuation(1.*cm),self.PhotoelectricFraction())

        return s 


if __name__ == '__main__':
    
    lyso = LYSO()

    print lyso  

    print "Average Lamda = %7.2f"%(lyso.AverageLamda()/nm)

    for z in drange(1., 11., 1.):
      print """
       z = %7.2g cm LYSO eff = %7.2g 
      """%(z,
           lyso.Efficiency(z*cm))
