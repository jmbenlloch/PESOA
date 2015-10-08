from LXe import *
from BLogging import *
import random as rnd
import numpy as np

log =logging.getLogger("LXSC")
log.setLevel(logging.INFO)
log.addHandler(ch)
deb = Switch(name="LXSC",switch=0)

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
               tDependenceGain=-0.008, 
               tDependenceDC=0.5,
               pulseRaise = 0.8*ns,
               pulseDecay= 1.3*ns,timeJitterFWHM=1*ns):

    """
    Defines a SiPM  :
    tDependenceDC = 0.5 for each 10 degree less 
    tDependenceGain 0.008 per degree less
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
    self.pulseRaise = pulseRaise
    self.pulseDecay = pulseDecay
    self.timeJitterFWHM = timeJitterFWHM

    self.activeBox = Box(activeSize,activeSize,1*mm)
    self.sensorBox = Box(sensorSize,sensorSize,1*mm)
  
  def Vendor(self):
    return self.vendor

  def Series(self):
    return self.series

  def Cost(self):
    return self.cost

  def SensorSize(self):
    return self.sensorSize

  def SensorArea(self):
    return self.sensorSize**2

  def ActiveArea(self):
    return self.activeSize**2

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

  def GainAtRoomT(self):
    return self.gain

  def DarkCurrentRateAtRoomT(self):
    return self.dCR

  def Gain(self,t):
    """
    Gain increases with temperature 
    """
    dt = celsius(STP_Temperature) + 20 - celsius(t)
    gf = dt*self.tDependenceGain
    return self.gain - gf*self.gain
  

  def DarkCurrent(self,t):
    """
    DCR decreases with temperature
    """
    dt = (celsius(STP_Temperature) + 20 - celsius(t))/10.
    gf = self.tDependenceDC**dt
    
    return self.DarkCurrentRateAtRoomT()*gf

  def PulseRaise(self):
    return self.pulseRaise

  def PulseDecay(self):
    return self.pulseDecay

  def TimeJitter(self):
    return self.timeJitterFWHM

  def __str__(self):
    s= """
      Vendor = %s 
      Series = %s
      Cost = %7.2f euro
      Sensor size = %7.2f mm  Active size = %7.2f mm 
      Microcell size = %7.2f micron 
      Number of microcells  = %7.2f PDE= %7.2f Gain= %7.4g cm
      Dark current rate (at 20 degrees) = %7.2f kHz
      Pulse Raise  = %7.2f ns PulseDecay = %7.2f ns TimeJitter= %7.4g ns
        """%(self.Vendor(), self.Series(), self.Cost()/euro,
          self.SensorSize()/mm,
          self.ActiveSize()/mm,
          self.MicroCellSize()/micron,
          self.MicroCellNumber(),
          self.PDE(),self.GainAtRoomT(),
          self.DarkCurrentRateAtRoomT()/kHz,
          self.PulseRaise()/ns,
          self.PulseDecay()/ns,
          self.TimeJitter()/ns)

    return s


class Sensor(SiPM):
  """
  Describes the sensor. It is a derived classs of SiPM adding the time response
  """


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
               tDependenceGain=-0.008, 
               tDependenceDC=0.5,
               pulseRaise = 0.8*ns,
               pulseDecay= 1.3*ns,timeJitterFWHM=0.1*ns,
               timeSeries =0, tauMax = 0):

    SiPM.__init__(self,series, vendor, cost, activeSize,sensorSize, 
               muCellSize,muCellNumber, muCellRecTime, pde, gain,
               dCR, CapacitanceAnodeToCathode, CapacitanceFastToCathode,
               tDependenceGain, tDependenceDC, pulseRaise,
               pulseDecay,timeJitterFWHM)

    self.STIME = timeSeries
    self.TAUMAX = tauMax

  def Jitter(self):
    """
    Time jitter: J is FWHM, divide by 2.355 to translate to sigma 
    """
    tj = rnd.gauss(0, self.TimeJitter()/2.355)
    return tj

  def SPE(self,t):
    """
    Single photoelectron response ( not normalized)
    """
    f = exp(-t/self.PulseRaise()) -exp(-t/self.PulseDecay())
    return f

  def SensorResponse(self,uvp):
    """
    Accumulated response to pe
    """

    # add the SiPM jitter to the time of the uvp
    t0 = uvp.Time() + self.Jitter()

    log.debug('+++Sensor Response +++')
    log.debug('t0 (ps)= %7.2f ', t0/ps)
    deb.Wait()

    if t0 > self.TAUMAX:  #loop out outliers
      return np.zeros(1)

    sper=[]
    nbins=0

    for tt in self.STIME:
      t= tt*ps #tt is expressed in ps
      if t > t0:
        spe = self.SPE(t)
      else:
        spe = 0.
      sper.append(spe)

      log.debug('t0 (ps)= %7.2f t (ps), spe', t0/ps, t/ps, spe)

      nbins+=1

    asper = np.array(sper)
    intg = sum(asper)
    asper = self.PDE()*asper/intg

    log.debug('Integral= %7.2f length of sper vector', intg, len(asper))
    log.debug('sper vector', asper)
    deb.Wait()

    return asper



class WLS:
  """
    Defines a Wavelength shifter
  """
  def __init__(self,name="TPB", lamda=420*nm, tau=1*ns):
    self.name = name
    self.lamda = lamda
    self.tau = tau

  def Name(self):
    return self.name

  def ScintillationWavelength(self):
    return self.lamda

  def Lifetime(self):
    return self.tau 

  def __str__(self):
    s= """
      Name = %s 
      Scintillation Wavelength of shifted light = %7.2f nm
      Lifetime for shifting = %7.2f ns
        """%(self.Name(), self.ScintillationWavelength()/nm, self.Lifetime()/ns)

    return s

class PLXSC:
  """
    Defines the Performance of the LXSC: 
    sigma_0 is the intrinsic resolution to add in cuadrature to the photoelectron statistics
    sensorEff is the detection efficiency of the sensor plane(s)
    uvRef is the reflectivity to the vuv light
    wlsRef is the reflectivity to the wls shifted (blue) light
    wlsEff is the WLS efficiency 
  """
  def __init__(self, sigmax=1*mm, sigmay=1*mm, sigmaz=1*mm, sigma0=0.05,
               sensorEff=0.9, wlsEff=0.8, uvRef=0.95, wlsRef=0.98):
    self.sigmax = sigmax 
    self.sigmay = sigmay
    self.sigmaz = sigmaz
    self.sigma0 = sigma0
    self.sensorEff = sensorEff
    self.wlsEff = wlsEff
    self.uvRef = uvRef
    self.wlsRef = wlsRef

  def SigmaX(self):
    return self.sigmax
  def SigmaY(self):
    return self.sigmay
  def SigmaZ(self):
    return self.sigmaz
  def Sigma0(self):
    return self.sigma0
  def SensorEfficiency(self):
    return self.sensorEff
  def WLSEfficiency(self):
    return self.wlsEff
  def ReflectivityUV(self):
    return self.uvRef
  def ReflectivityWLS(self):
    return self.wlsRef


  def __str__(self):
    s= """
      sigmax=%7.2f mm, sigmay=%7.2f mm, sigmaz=%7.2f ,sigma0=%7.2f
      sensorEff=%7.2f, wlsEff=%7.2f, uvRef=%7.2f, wlsEff=%7.2f
      
        """%(self.SigmaX()/mm, self.SigmaY()/mm, self.SigmaZ()/mm, self.Sigma0(),
             self.SensorEfficiency(), self.WLSEfficiency(), self.ReflectivityUV(), 
             self.ReflectivityWLS())

    return s
  
class LXSC:
    """
    Defines a Liquid Xenon Scintillating Cell 
    instMask is a list which defines the mask of instrumented faces.
    it runs from 0 to 2: 
      0 = entry face/exit face (x), 
      1=  left face/right face (y) 
      2=  topface/bottom face (z), 
    Faces are assumed to be instrumented in pairs 
    [1,0,0] means that both x faces are instrumented and the other faces are not
    """
    def __init__(self,lxe, box, wls, plxsc, sensor, pitch,
                 instMask=[1,0,0]):

      self.lxe = lxe
      self.box=box
      self.plxsc = plxsc
      self.sipm=sensor
      self.wls = wls
      self.pitch = pitch 
      self.instMask = instMask
      
    def LXe(self):
      return self.lxe

    def Box(self):
      return self.box
    
    def WLS(self):
      return self.wls

    def PLXSC(self):
      return self.plxsc

    def SensorPitch(self):
      return self.pitch

    def Sensor(self):
      return self.sipm

    def Mass(self):
      return self.Box().V()*self.LXe().Density() 

    def CostOfXenon(self):
      return float((self.Mass()/g)*self.LXe().CostPerGram())

    def ScintillationPhotons(self,E):
      return self.LXe.ScintillationPhotons(E)

    def SPhotonsAt511KeV(self,i):
      return self.LXe.SPhotonsAt511KeV(i)

    def InstrumentedMask(self):
      return self.instMask

    def SiXZ(self):
      return int(self.Box().XZ()/self.SensorPitch()**2)-1

    def SiYZ(self):
      return int(self.Box().YZ()/self.SensorPitch()**2)-1

    def SiXY(self):
      return int(self.Box().XY()/self.SensorPitch()**2)-1

    def NumberOfSiPM(self):
      nsipm = 2*self.SiXY()*self.instMask[0] 
      nsipm += 2*self.SiYZ()*self.instMask[1]
      nsipm +=2*self.SiXZ()*self.instMask[2] 

      return int(nsipm)

    def NX(self):
      n = self.Box().X()/self.SensorPitch()
      return (int)(n)

    def NY(self):
      n = self.Box().Y()/self.SensorPitch()
      return (int)(n)
  
    def Nz(self):
      n = self.Box().Z()/self.SensorPitch()
      return (int)(n)

    def CostOfSiPM(self):
      return float(self.NumberOfSiPM()*self.Sensor().Cost())

    def CostOfCell(self):
      return float(self.CostOfSiPM()+self.CostOfXenon())

    def SensorCoverage(self):
      return self.Sensor().ActiveArea()/self.SensorPitch()**2
    
    
    def __str__(self):
        
        s= """
        LXe =%s
        BOX = %s
        WLS = %s
        Performance LXSC  = %s
        SiPM = %s
        mask = %s
        pitch = %7.2f mm
        mass = %7.2f g
        number of SiPMs = %d
        cost of Lxe = %7.2f euro
        cost of SiPM = %7.2f euro
        cost of Cell = %7.2f euro

        
       
    """%(self.LXe(),self.Box(),self.WLS(),self.PLXSC(),self.Sensor(),
         self.InstrumentedMask(),
         self.SensorPitch(), 
         self.Mass()/g,
         self.NumberOfSiPM(),
         self.CostOfXenon(), 
         self.CostOfSiPM(),
         self.CostOfCell())
         

        return s

if __name__ == '__main__':
    
    lxe = LXe()
    box = Box(25*mm,25*mm,30*mm)
    tpb = WLS()
    plxsc = PLXSC()
    sipm = SiPM()
    lxsc = LXSC(lxe,box,tpb,plxsc,sipm,6.2*mm,[1,0,0])
    
    print lxsc
    
    print """dark current of SiPM:
              at room temperature = %7.2f kHz 
              at LXe temperature (%7.2f Kelvin) = %7.2f kHz

    """%(
      sipm.DarkCurrentRateAtRoomT()/kHz,
      lxe.TemperatureAt1Bar()/kelvin,
      sipm.DarkCurrent(lxe.TemperatureAt1Bar())/kHz)

    print """Gain of SiPM:
              at room temperature = %7.2g 
              at LXe temperature (%7.2f Kelvin) = %7.2g
              
    """%(
      sipm.GainAtRoomT(),
      lxe.TemperatureAt1Bar()/kelvin,
      sipm.Gain(lxe.TemperatureAt1Bar()))

    print """
      LXSC efficiency in sensor plane = %7.2f
    """%lxsc.SensorCoverage()


