from Util import *
from LXe import *

class Box:

  def __init__(self,x,y,z):
    """
    Defines a box of dimensions x,y,z
    z is considered the longitudinal dimension 
    """
    self.x=x
    self.y=y
    self.z=z
        
  def V(self):
    return self.x*self.y*self.z

  def X(self):
    return self.x

  def Y(self):
    return self.x

  def Z(self):
    return self.z

  def XZ(self):
    return self.x*self.z

  def XY(self):
    return self.x*self.y

  def YZ(self):
    return self.z*self.y
    

  def __str__(self):
        
    s= """
        Box:
        x = %7.2g cm y = %7.2g cm z = %7.2g cm
        sxy = %7.2g cm2 syz = %7.2g cm2 sxz = %7.2g cm2
        V = %7.2g cm3
      """%(self.X()/cm, self.Y()/cm, self.Z()/cm,
           self.XY()/cm2, self.YZ()/cm2, self.YZ()/cm2, 
           self.V()/cm3)
    return s

class SiPM:

  def __init__(self,series="J", vendor="SENSL",
               cost=20*euro,
               activeSize=6*mm,sensorSize=6.13*mm, 
               muCellSize=35*micron,
               muCellNumber=22300,
               muCellRecTime=210*ns,
               pde=0.51,
               gain=6e+6,
               dCR=70*mm2*kHz,
               CapacitanceAnodeToCathode=3400*pF,
               CapacitanceFastToCathode=48*pF,
               tDependenceGain=-0.008/kelvin,
               tDependenceDC=0.05/kelvin):
    """
    Defines a SiPM  
    """

    self.series=series
    self.vendor=vendor
    self.cost = cost
    self.activeSize=activeSize
    self.sensorSize=sensorSize
    self.muCellSize=muCellSize
    self.muCellNumber=muCellNumber
    self.muCellRecTime=muCellRecTime
    self.pde=pde
    self.gain=gain
    self.dCR=dCR*(self.activeSize/mm)**2
    self.CapacitanceAnodeToCathode=CapacitanceAnodeToCathode
    self.CapacitanceFastToCathode=CapacitanceFastToCathode
    self.tDependenceGain=tDependenceGain
    self.tDependenceDC=tDependenceDC

    self.activeBox = Box(activeSize,activeSize,activeSize)
    self.sensorBox = Box(sensorSize,sensorSize,sensorSize)
  
  def Vendor(self):
    return self.vendor

  def Series(self):
    return self.series

  def Cost(self):
    return self.cost

  def SensorSize(self):
    return self.sensorSize

  def ActiveSize(self):
    return self.activeSize

  def SensorBox(self):
    return self.sensorBox

  def ActiveBox(self):
    return self.activeBox

  def MicroCellSize(self):
    return self.muCellSize

  def MicroCellNumber(self):
    return self.muCellNumber

  def MicroCellRecoveryTime(self):
    return self.muCellRecTime

  def PDE(self):
    return self.pde

  def Gain(self):
    return self.gain

  def DarkCurrentRate(self):
    return self.dCR

  def __str__(self):
    s= """
      Vendor = %s 
      Series = %s
      Cost = %7.2f euro
      Sensor size = %7.2f mm  Active size = %7.2f mm 
      Microcell size = %7.2f micron 
      Number of microcells  = %7.2f PDE= %7.2f Gain= %7.4g cm
      Dark current rate (at 20 degrees) = %7.2f kHz
        """%(self.Vendor(), self.Series(), self.Cost()/euro,
          self.SensorSize()/mm,
          self.ActiveSize()/mm,
          self.MicroCellSize()/micron,
          self.MicroCellNumber(),
          self.PDE(),self.Gain(),
          self.DarkCurrentRate()/kHz)

    return s


class LXSC:
    """
    Defines a Liquid Xenon Scintillating Cell 
    instMask is a list which defines the mask o instrumented faces.
    it runs from 0 to 5: (0 = entry face, 1= exit face, 2= left face, 
    3 = right face, 4 = upper face, 5 = bottom face): 1 means face
    is instrumented with SiPMs, 0 means is not.
    """
    def __init__(self,Xe, box, siPM, pitch,
                 instMask=[1,1,0,0,0,0]):

        self.Xe = Xe
        self.box=box
        self.siPM=siPM
        self.pitch = pitch 
        self.instMask = instMask

    def SensorPitch(self):
      return self.pitch

    def Box(self):
      return self.box

    def SiPM(self):
      return self.siPM

    def Mass(self):
        return self.box.V()*self.Xe.Density() 

    def CostOfXenon(self):
      return float((self.Mass()/g)*self.Xe.CostPerGram())

    def ScintillationPhotons(self,E):
        return self.Xe.ScintillationPhotons(E)

    def InstrumentedMask(self):
        return self.instMask

    def SiXZ(self):
        return int(self.Box().XZ()/self.SensorPitch()**2)-1

    def SiYZ(self):
        return int(self.Box().YZ()/self.SensorPitch()**2)-1

    def SiXY(self):
        return int(self.Box().XY()/self.SensorPitch()**2)-1

    def NumberOfSiPM(self):
      nsipm = self.SiXY()*self.instMask[0] + self.SiXY()*self.instMask[1]
      nsipm += self.SiYZ()*self.instMask[2]+ self.SiYZ()*self.instMask[3]
      nsipm +=self.SiXZ()*self.instMask[4] +self.SiXZ()*self.instMask[5]

      return int(nsipm)

    def CostOfSiPM(self):
      return float(self.NumberOfSiPM()*self.SiPM().Cost())

    def CostOfCell(self):
      return float(self.CostOfSiPM()+self.CostOfXenon())
     

    def __str__(self):
        
        s= """
        LXSC
        BOX = %s
        SiPM = %s
        pitch = %7.2f mm
        mass = %7.2f g
        number of SiPMs = %d
        cost of Lxe = %7.2f euro
        cost of SiPM = %7.2f euro
        cost of Cell = %7.2f euro

        
       
    """%(self.Box(),
         self.SiPM(),
         self.SensorPitch(), 
         self.Mass()/g,
         self.NumberOfSiPM(),
         self.CostOfXenon(), 
         self.CostOfSiPM(),
         self.CostOfCell())
         

        return s

if __name__ == '__main__':
    
    lxe = LXe()
    box = Box(50*mm,50*mm,50*mm)
    sipm = SiPM()
    lxsc = LXSC(lxe,box,sipm,6.2*mm,[1,1,0,0,0,0])
    
    print lxsc