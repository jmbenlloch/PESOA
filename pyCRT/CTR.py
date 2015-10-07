from Photon import *

#import pylab as P

TAUMAX = 3*TAU1
TBINS = 0.01*ns
NBINS = (int) (TAUMAX/TBINS)

TLEVEL = 1 # in pes

MAX_BOUNCES = 100




if __name__ == '__main__':

	c1 = TCanvas( 'c1', 'CTR', 200, 10, 600, 800 )
	

	wbox = Box(LX,LY,LZ)
	pt = PhotonTransport(wbox)
	pg = PhotonGenerator(wbox)
	sipm =Sensor(PDE=SIPM_PDE, J=SIPM_J, taur=SIPM_PR,taud=SIPM_PD)
	sce = ScintillationEvent(wbox,pg,sipm,TAU1,NPHOTONS,TLEVEL)

	ht0 = TH1F("ht0", "photon creation time (ns)", NBINS, 0., 3*TAUMAX/ns)
	htime = TH1F("htime", "photon time (ns)", NBINS, 0., 3*TAUMAX/ns)
	hpath = TH1F("hpath", "photon path (mm)", 100, 0., 500.)
	hb = TH1F("hb", "number of bounces", 20, 0., 20.)
	hspe = TH1F("hspe", "spe(au)", NBINS, 0.,3*TAUMAX/ns)
	hCTR1 = TH1F("hCTR1", "ctr1", 200, 800.,1000.)
	hCTR2 = TH1F("hCTR2", "ctr2", 50, 1800.,2200.)
	

	print """
		CTR calculation
		number of events = %d 
		using tau = %7.2f ns
		taumax = %7.2f ns
		number of photons for this tau = %d
		number of time bins for photons = %d
		number of time bins for SPE = %d
		Box dimensions: x = %7.2f mm y = %7.2f mm z = %7.2f mm 
		SiPM PDE = %7.2f
		SiPM Jitter = %7.2f
		SiPM raise constant %7.2f ns
		SiPM decay constant %7.2f ns
		Trigger level = %7.2f
			
	"""%(NEVENTS,TAU/ns,TAUMAX/ns,NG_tau1,NBINS,2*NBINS,
			LX/mm,LY/mm,LZ/mm,SIPM_PDE,SIPM_J,SIPM_PR/ns,SIPM_PD/ns,TLEVEL)

	wait()


	for iev in range(NEVENTS):	
		print "event ", iev
		
		sce.Event()
		sce.Trigger()
		
	c1.Divide(2,2)
	c1.cd(1)
	gPad.SetLogy()
	ht0.Draw()
	
	c1.cd(2)
	gPad.SetLogy()
	htime.Draw()

	c1.cd(3)
	gPad.SetLogy()
	hpath.Draw()

	c1.cd(4)
	gPad.SetLogy(0)
	hCTR2.Draw()

	c1.Show()
	wait()

	# c1.Divide(1,2)
	# c1.cd(1)
	# gPad.SetLogy(0)
	# hCTR1.Draw()
	
	# c1.cd(2)
	# gPad.SetLogy(0)
	# hCTR2.Draw()

	# c1.Show()
	# wait()
	

	#FPULSE.Draw()
	
	# FTAU3.Draw()
	# FTAU2.Draw("same")
	# FTAU1.Draw("same")
	# c1.Show()
	# wait()
	

	
	


