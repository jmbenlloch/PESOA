from Event import *
from THisto import *
import random as rnd

log =logging.getLogger("CTR")
log.setLevel(logging.INFO)
log.addHandler(ch)
deb = Switch(name="CTRDebug",switch=0)
hsw = Switch(name="CTRHistos",switch=1)


def Histos(th,lxsc,tauIndex):
	lxe = lxsc.LXe()
	wbox =lxsc.Box()
	nx = lxsc.NX()
	ny = lxsc.NY()

	if hsw.Switch() == 1:
		th.BookH1("ht0", "photon t0 (ns)", 200, 0., 5*lxe.Lifetime(tauIndex)/ns)
		th.BookH1("htime", "photon time (ns)", 200, 0., 5*lxe.Lifetime(tauIndex)/ns)
		th.BookH1("hctime", "photon corrected time true (ps)", 200, 0., 1000)
		th.BookH1("hcrtime", "photon corrected time rec (ps)", 200, 0., 1000)

		th.BookH1("hctimeR", "time - ctime (ns)", 200, 0., 5*lxe.Lifetime(tauIndex)/ns)
		th.BookH1("hcrtimeR", "time - crtime (ns)", 200, 0., 5*lxe.Lifetime(tauIndex)/ns)
	
		th.BookH1("hx0", "x0 (mm) ", 50, 0., wbox.X()/mm)
		th.BookH1("hy0", "y0 (mm) ", 50, 0., wbox.Y()/mm)
		th.BookH1("hz0", "z0 (mm) ", 50, 0., wbox.Z()/mm)

		th.BookH1("hz", "z (mm) ", 50, 0., (wbox.Z()+10)/mm)
		th.BookH1("hdph", " distance photon true (mm) ", 50, 0., wbox.Z()/mm)
		th.BookH1("hdphr", " distance photon rec (mm) ", 50, 0., wbox.Z()/mm)		
							
		th.BookH1("ht0fe", "t0 first pe (ps)", 200, 0., 100)
		th.BookH1("htfe", "time first pe (ps)", 200, 0., 500)
		th.BookH1("hctfe", "ctime first pe (ps)", 200, 0., 500)
		th.BookH1("hcrtfe", "crtime first pe (ps)", 200, 0., 500)

		th.BookH1("hfctimeR", "timeFirstPe - ctimeFirstPe (ps)", 200, -200., 200)
		th.BookH1("hfcrtimeR", "timeFirstPe - crtimeFirstPe", 200, -200., 200)

		th.BookH1("hctimeRR", "ctime - crtime", 100, -100., 100)
		th.BookH1("hfctimeRR", "first: ctime - crtime", 100, -100., 100)
		th.BookH1("hxpes", "x pes in mm", nx, 0., wbox.X()/mm)
		th.BookH1("hypes", "y pes in mm", ny, 0., wbox.Y()/mm)
		th.BookH2("hxypes", "x-y pes in mm", nx, 0., wbox.X()/mm,ny, 0., wbox.Y()/mm)


def EventVUV(tauIndex,PDE,NEVENTS):
	"""
	Generates an event 
	"""

	lxe = LXe()
	sipm = Sensor(timeSeries =0, tauMax = 0,
		          pde=PDE,pulseRaise = 1.0*ns, pulseDecay= 3.2*ns,
				  timeJitterFWHM=0*ns)

	sigma0 = Point3D(x=1*mm,y=1*mm,z=1*mm) # resolution in interaction
	sigma = Point3D(x=2*mm,y=2*mm,z=1*mm)  # resolution in SiPM
	SensorEff=0.9
	WlsEff=0.8 
	UvRef=0.95 
	WlsRef=0.98

	plxsc = PLXSC(sigmax=sigma0.x, sigmay=sigma0.y, sigmaz=sigma0.z, sigma0=0.05,
                  sensorEff=SensorEff, wlsEff=1, uvRef=UvRef, wlsRef=1)

	wbox = Box(50*mm,50*mm,50*mm)
	tpb = WLS()

	imask = [0,0,1]
	lxsc = LXSC(lxe,wbox,tpb,plxsc,sipm,6.2*mm,imask)

	histoPath = "../histos/CTRHistos_VUV_tau_%d_PDE_%7.2f_BOX_%7.2f"%(tauIndex,PDE,wbox.X()/cm)
	fHistoFile = TFile(histoPath,"RECREATE")

	 
	th = THisto()
	Histos(th,lxsc,tauIndex)
	
	sce = ScintillationEvent(lxsc,tauIndex)
	
	print """
		CTR calculation: VUV photons
		number of events = %d 
		using tau = %7.2f ns
		number of photons for this tau = %d
		lxsc = %s
		
			
	"""%(NEVENTS,lxe.Lifetime(tauIndex)/ns,sce.np,lxsc)

	wait()

	for iev in range(NEVENTS):	
		x = rnd.uniform(0.,wbox.X())
		y = rnd.uniform(0.,wbox.Y())
		z = rnd.uniform(0.,wbox.Z())
		w = lxe.Efficiency(511*keV,z)
		
		log.info('event  = %d: z = %7.2f, w = %7.2f z*w %7.2f',iev, z/mm,w/mm,z*w/mm)
		if hsw.Switch() == 1:
			th.FillH1("hx0",x/mm)
			th.FillH1("hy0",y/mm)
			th.FillH1("hz0",z*w/mm)
		#wait()

		#ixyz = Point3D(wbox.X()/2,wbox.Y()/2,wbox.Z()/2)
		ixyz = Point3D(x,y,z*w)
		txyz = Point3D(-9999,-9999,-9999)

		log.debug('ixyz  = %s: txyz = %s',ixyz, txyz)

		PES = sce.Event(ixyz,txyz)
		
		ipes = 0
		t0FirstPe = 1*second
		timeFirstPe = 1*second
		ctimeFirstPe = 1*second
		crtimeFirstPe = 1*second

		w = 1./len(PES)

		for pes in PES:
			ipes+=1
			log.debug('pes number  = %d pes = %s ',
				iev,pes)
			photon = pes.photon
			x0 = photon.xyzt0.x
			y0 = photon.xyzt0.y
			z0 = photon.xyzt0.z
			t0 = photon.xyzt0.t

			log.debug('photon init coor: x0 = %7.2f mm, y0 = %7.2f mm, z0= %7.2f  mm, t0= %7.2f ns ',
				 x0/mm,y0/mm,z0/mm,t0/ns)

			x = photon.xyzt.x
			y = photon.xyzt.y
			z = photon.xyzt.z
			t = photon.xyzt.t 

			log.debug('photon final (pes) coor: x = %7.2f mm, y = %7.2f mm, z= %7.2f mm, t0= %7.2f ns ',
				 x/mm,y/mm,z/mm,t/ns)

			if hsw.Switch() == 1:
				th.FillH1("hxpes", x/mm)
				th.FillH1("hypes", y/mm)
				th.FillH2("hxypes",x/mm,y/mm)
				th.FillH1("ht0", t0/ns)
				th.FillH1("htime",t/ns)
			
			
			#Compute corrected time
			# we measure time stamp of pes t.
			# ctime is an estimator of t0 which is computed as tc = d * n/c
			# where the distance d:
			# d = sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)
			# (x0,y0,z0) are the interaccion coordiantes (measured with some error)
			# (x,y,z) are the detection coordinates (estimated from SiPM position) 

			x0r = x0 + rnd.gauss(0, sigma0.x) # reconstructed x0
			y0r = y0 + rnd.gauss(0, sigma0.y) # reconstructed y0
			z0r = z0 + rnd.gauss(0, sigma0.z) # reconstructed z0

			xr = x + rnd.gauss(0, sigma.x) # reconstructed x
			yr = y + rnd.gauss(0, sigma.y) # reconstructed y
			zr = z + rnd.gauss(0, sigma.z) # reconstructed z

			dr = sqrt((xr-x0r)**2 + (yr-y0r)**2 + (zr-z0r)**2)
			d = sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
			
			ctime = d*photon.n/c_light  # true corrected time
			crtime = dr*photon.n/c_light  # reconstructed corrected time

			pes.ctime = ctime
			pes.crtime = crtime 

			log.debug('dr = %7.2f mm, d = %7.2f mm',
				 dr/mm,d/mm)
			log.debug('t0 = %7.2f ns, t = %7.2f ns, ctime = %7.2f ns, crtime = %7.2f ns',
				 pes.photon.T0()/ns,pes.photon.Time()/ns,ctime/ns, crtime/ns)

			log.debug('time - ctime  = %7.2f ns, time -crtime = %7.2f ns, ctime - crtime %7.2f ps',
			(pes.photon.T0() - ctime)/ns,(pes.photon.Time() - crtime)/ns, (ctime - crtime)/ps)
			
			if hsw.Switch() == 1:
				th.FillH1("hdph", d/mm)
				th.FillH1("hdphr", dr/mm)
				th.FillH1("hz",z/mm)

				th.FillH1("hctime", ctime/ps)
				th.FillH1("hcrtime", crtime/ps)

				th.FillH1("hctimeR", (pes.photon.Time() - ctime)/ns)
				th.FillH1("hcrtimeR", (pes.photon.Time() - crtime)/ns)
				th.FillH1("hctimeRR", (ctime - crtime)/ps)
				
			pesTime = pes.photon.Time() - ctime 			

			if pes.photon.Time() < timeFirstPe:
				timeFirstPe = pes.photon.Time()
				t0FirstPe = pes.photon.T0()
				ctimeFirstPe = ctime
				crtimeFirstPe = crtime
			deb.Wait()

		if hsw.Switch() == 1:
			th.FillH1("ht0fe", t0FirstPe/ps)
			th.FillH1("htfe", timeFirstPe/ps)
			th.FillH1("hctfe", ctimeFirstPe/ps)
			th.FillH1("hcrtfe", crtimeFirstPe/ps)
			th.FillH1("hfctimeR", (timeFirstPe - ctimeFirstPe)/ps)
			th.FillH1("hfcrtimeR", (timeFirstPe - crtimeFirstPe)/ps)
			th.FillH1("hfctimeRR", (ctimeFirstPe - crtimeFirstPe)/ps)
			
		log.info('Number of PES = %d',len(PES))
		log.info('first pes: t0 = %7.2f ps, time = %7.2f ps, ctime = %7.2f, ps crtime = %7.2f ps  ',
			t0FirstPe/ps,timeFirstPe/ps,ctimeFirstPe/ps,crtimeFirstPe/ps)
		log.info('time - ctime  = %7.2f ps, time -crtime = %7.2f ps, ctime - crtime %7.2f ps',
			(timeFirstPe - ctimeFirstPe)/ps,(timeFirstPe - crtimeFirstPe)/ps, 
			(ctimeFirstPe - crtimeFirstPe)/ps)

		#wait()
		

		#sce.Trigger()
	if hsw.Switch() == 1:	
		th.DrawList(["ht0","htime","hctime","hcrtime"],xd=2,yd=2,yscale="log")
		th.DrawList(["hctimeR","hcrtimeR"],xd=1,yd=2,yscale="log")
		th.DrawList(["hx0","hy0","hz0"],xd=1,yd=3,yscale="lin")
		th.DrawList(["hz","hdph","hdphr"],xd=1,yd=3,yscale="lin")
		th.DrawList(["hxpes","hypes","hxypes"],xd=1,yd=3,yscale="lin")
		th.DrawList(["ht0fe","htfe","hctfe","hcrtfe"],xd=2,yd=2,yscale="lin")
		th.DrawList(["hfctimeR","hfcrtimeR","hctimeRR","hfctimeRR"],xd=2,yd=2,yscale="lin")
	
	
	
if __name__ == '__main__':
	NEVENTS = 1000
	
	EventVUV(1,0.2,NEVENTS)
	
