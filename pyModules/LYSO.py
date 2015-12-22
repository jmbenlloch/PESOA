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

    for z in drange(1., 11., 1.):
      print """
       z = %7.2g cm LYSO eff = %7.2g 
      """%(z,
           lyso.Efficiency(z*cm))
