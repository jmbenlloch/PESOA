from abc import ABCMeta, abstractmethod
from Util import *

# class Cell:
#   __metaclass__ = ABCMeta
#   @abstractmethod
#   def Volume(self):
#     pass

#   @abstractmethod
#   def Sr(self):
#     pass

#   @abstractmethod
#   def L(self):
#     pass

#   @abstractmethod
#   def R(self):
#     pass

class XYZBox:

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

  def Sxy(self):
    return self.x*self.y

  def Sxz(self):
    return self.x*self.z

  def Syz(self):
    return self.x*self.z


  def __str__(self):
        
    s= """
        XYZBox:
        x = %7.2f y = %7.2f  z = %7.2f
        Sxy = %7.2f Sxz = %7.2f Syx = %7.2f 
        V = %7.2f
      """%(self.x, self.y, self.z, self.Sxy(),
        self.Sxz(), self.Syz(), self.V())
    return s

class Box(XYZBox):

  def __init__(self, boxCoord):
    """
        Defines a box in terms of the box coordinates in the reference system
        of G4. The box coordinates are defined by 8 vertices.
        Each vertex is an (x,y,z) point. 
        Define the box with the following convention:
        v1 = (xl,yb,zb), where xl = lefmost x coordinate,
                               yb = bottom y coordinate
                               zb = back z coordinate
        x-------x (xr,yt)
        |       |
        |       |
        x-------x            zb ----- zf
        (xl,yb)

        v2 = (xl,yt,zb), where xl = lefmost x coordinate,
                               yt = top y coordinate
                          

        v3 = (xr,yb,zb), where xr = rightmost x coordinate
        v4 = (xr,yt,zb)
        v5 = (xl,yb,zf), where zf = front z coordinate
        v6 = (xl,yt,zf)
        v7 = (xl,yt,zf)
        v8 = (xr,yb,zf)                           
 
    """
    self.boxCoord =boxCoord
    v1 = self.boxCoord[0] 
      
    self.xmin = v1[0]
    self.ymin = v1[1]
    self.zmin =v1[2]

    v8 = self.boxCoord[7] 
    
    self.xmax = v8[0]
    self.ymax = v8[1]
    self.zmax = v8[2]

    x = self.xmax-self.xmin
    y = self.ymax-self.ymin
    z = self.zmax-self.zmin

    XYZBox.__init__(self, x,y,z)

  def Active(self,coord):
    """
    returns true if coord = (x,y,z) is in the active volume of the box
    defined by xmin-xmax, ymin-ymax, zmin-zmax, false otherwise
    """

    x,y,z = coord
    box = False
    if x > self.xmin and x < self.xmax:
      if y > self.xmin and y < self.ymax:
        if z > self.zmin and z < self.zmax:
          box = True
    
    return box
         

  def __str__(self):
       
    s =  XYZBox.__str__(self)
    s+= """
        Box:
        xmin = %7.2f ymin = %7.2f zmin = %7.2f
        xmax = %7.2f ymax = %7.2f zmax = %7.2f
      """%(self.xmin, self.ymin, self.zmin,
           self.xmax, self.ymax, self.zmax)
    return s

if __name__ == '__main__':
    
    boxCoord =[[-12.8, -12.8, -130],[-12.8, 12.8, -130],[12.8, -12.8, -130],  
    [12.8, 12.8, -130],[-12.8, -12.8, -100],[-12.8, 12.8, -100],
    [12.8, -12.8, -100],[12.8, 12.8, -100]]
    box = Box(boxCoord)
    fid = box.Active([10.,5., -110.])
    print fid
    
    print box