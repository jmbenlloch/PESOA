from math import *

D=220 #mm
n = 22 #number of modules in 2pi
nc = n/4. #number of modules per cuadrant
theta = 2*pi/n  #segmentation angle
L = 30 #mm depths of scanner

D=800 #mm
n = 50 #number of modules in 2pi
nc = n/4. #number of modules per cuadrant
theta = 2*pi/n  #segmentation angle
L = 20 #mm depths of scanner

display = False


#def Cost():



def ModuleVolume(mi,mo,L):
	"""
	The module is a pyramidal fustrum (a truncated pyramid) of small area
	A1 = mi^2 and large area A2 = mo^2 the height is L
	the volume is 1/3h(A1+A2+Sqrt(A1A2))
	"""
	return (1./3.)*L*(mi*mi + mo*mo + sqrt(mi*mo))

def AppendP(P,x,y):
	P.append((x,y))

def PrintP(P):
	print "(X,Y) in mm"
	for p in P:
		x = p[0]
		y = p[1]
		print """
			(%7.2f,%7.2f)
			"""%(x,y)

def FillP(theta,R,P):
	for n in range(nc+1):
		t = n*theta
		x = R*cos(t)
		y = R*sin(t)
		AppendP(P,x,y)
		
	for n in range(nc+1):
		p =P[n]
		if p[0] != 0: 
			x = -p[0]
			y = p[1]
			AppendP(P,x,y)

	for n in range(nc):
		p =P[n]
		if p[1] != 0:
			x = p[0]
			y = -p[1]
			AppendP(P,x,y)

	for n in range(nc):
		p=P[n]
		if p[0] != 0 and p[1] != 0:
			x = -p[0]
			y = -p[1]
			AppendP(P,x,y)

def Display(P,RI,R):
	print """
			cI = Graphics[Circle[{0, 0}, %7.2F]]
			"""%(RI)
	print """
			cO = Graphics[Circle[{0, 0}, %7.2F]]
			"""%(R)

	n=1
	for p in P:
		print """
			l%d = Graphics[Line[{{0, 0}, {%7.2F, %7.2F}}]]
			"""%(n,p[0],p[1])
		n+=1

	lines=[]
	for i in range(len(P)):
		n=i+1
		lines.append("l"+str(n))
	
	for l in lines:
		print l,

if __name__ == '__main__':
	RI=D/2.
	RO = RI + L

	cI = 2*pi*RI
	cO = 2*pi*RO

	mI = cI/n
	mO = cO/n

	V = ModuleVolume(mI/10.,mO/10.,L/10.)
	print """
	PET internal diameter = %7.2f mm
	PET internal radius = %7.2f mm
	PET external radius = %7.2f mm
	number of modules = %d
	theta = %7.2f radians

	mI for module = %7.2f mm
	mO for module = %7.2f mm
	Volume of module %7.2f cm3
	"""%(D,RI,RO,n,theta,mI,mO,V)

	if display == True:
		P=[]
		FillP(theta,R,P)
		Display(P,RI,R)


