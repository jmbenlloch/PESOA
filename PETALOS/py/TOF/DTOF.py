from Centella.AAlgo import AAlgo

from Util import *
from Geometry import *
from TOF import *
from Scintillator import *

"""
This algorithm computes the Coincidence Resolution Time. 
"""


class DTOF(AAlgo):

	############################################################
	def __init__(self, param=False, level=1, label="", **kargs):

		"""
		DTOF Algorithm
		"""
		#self.m.log(1, 'Constructor()')

		### GENERAL STUFF
		self.name = 'DTOF'
		
		AAlgo.__init__(self, param, level, self.name, 0, label, kargs)

    ### PARAMETERS
		
		
		self.debug = self.ints["Debug"]  #used to stop program at key break points
		self.SCINT = self.strings["SCINTILLATOR"]

		#load coordinates of box and fiducial box

		if self.SCINT == "LXE":
			self.scint = LXe() #lxe properties
		elif self.SCINT == "LYSO":
			self.scint = LYSO()
		else:
			print "scintillator not yet implemented"
			sys.exit()

		print self.scint


		if self.debug == 1:
			wait()


	############################################################		
	def initialize(self):

		self.m.log(1, 'Initialize()')
		
		### Defining histos
		# Event energy histogram

		self.Histos()

    ### Counters:
		self.numInputEvents = 0
		self.numOutputEvents = 0

		
		return


	############################################################
	def execute(self, event=""):

		self.m.log(2, 'Execute()')	

		#get time map from previous algo

		self.timeMap = self.logman["USER"].gparam["TimeMap"]
		
		#Compute DT
		self.ComputeDTOF() 
		
		
		if self.debug == 1:
			wait()	

		
		return True

############################################################
	def finalize(self):

		self.m.log(1, 'Finalize()')

		self.m.log(1, 'Input  Events: ', self.numInputEvents)
		self.m.log(1, 'Output Events: ', self.numOutputEvents)
		

		self.logman["USER"].ints[self.alabel("InputEvents")] = self.numInputEvents
		self.logman["USER"].ints[self.alabel("OutputEvents")] = self.numOutputEvents
	
		#self.tman.save(file_name=self.FTREES)

		return


###########################################################
	def ComputeDTOF(self):
		"""
		Compute DT for the first PE using directly gamma information
		DTOF = (1/2)*(dt - dtg - dpg)
		where: DT is the difference in TOF
		dt = t2 -t1 difference of time stamp of first PE 
		dtg = tg2 -tg1 difference of time stamp of gamma 1 and tGamma2 
		dpg = difference of path1 and path2, where 
		path = (n/c)*d and distance is the distance between gamma and SiPM vertex. 
		"""

		self.m.log(2," ---ComputeDTOF---")

		vertexBox1 =self.timeMap.InteractionVertex(box=1)
		vertexBox2 =self.timeMap.InteractionVertex(box=2)
		siPMHit1 = self.timeMap.SiPmHit(box=1,index=0)
		siPMHit2 = self.timeMap.SiPmHit(box=2,index=0)


		dt  = siPMHit2.TimeFirstPE() - siPMHit1.TimeFirstPE()
	
		dtg = vertexBox2.t - vertexBox1.t
		
		dbox1 = distance(siPMHit1.XYZ(),vertexBox1.XYZ())
		tpath1 = dbox1*self.scint.RefractionIndex()/c_light
		dbox2 = distance(siPMHit2.XYZ(),vertexBox2.XYZ())
		tpath2 = dbox2*self.scint.RefractionIndex()/c_light

		dpg = tpath2 - tpath1

		t1stPeBox1 = siPMHit1.TimeFirstPE() - tpath1 - vertexBox1.t
		t1stPeBox2 = siPMHit2.TimeFirstPE() - tpath2 - vertexBox2.t

		self.hman.fill(self.TimeFirstPeBox1_histo_name,t1stPeBox1/ps)
		self.hman.fill(self.TimeScintLightBox1_histo_name,tpath1/ps)
		self.hman.fill(self.TimeScintLightBox2_histo_name,tpath2/ps)
		self.hman.fill(self.TimeFirstPeBox2_histo_name,t1stPeBox2/ps)
		self.hman.fill(self.TimeFirstPeBox1Box2_histo_name,t1stPeBox1/ps,t1stPeBox2/ps)
			
		dtof = 0.5*(dt - dtg - dpg)
		self.hman.fill(self.dpg_histo_name,dpg/ps)
		self.hman.fill(self.dt_histo_name,dt/ps)
		self.hman.fill(self.dtg_histo_name,dtg/ps)
		self.hman.fill(self.dtof_histo_name,dtof/ps)
		self.hman.fill(self.dtof2_histo_name,dtof/ps)
		self.hman.fill(self.dtof3_histo_name,dtof/ps)
		
		
		self.m.log(2,
			"Box1: t1stPE =%7.2f ps, tgamma= %7.2f ps, tpath =%7.2f ps, dt1stPE =%7.2f ps "%(
				siPMHit1.TimeFirstPE()/ps,vertexBox1.t/ps,tpath1/ps, 
				(siPMHit1.TimeFirstPE()-vertexBox1.t-tpath1)/ps))
		self.m.log(2,
			"Box2: t1stPE =%7.2f ps, tgamma= %7.2f ps, tpath =%7.2f ps, dt1stPE =%7.2f ps "%(
				siPMHit2.TimeFirstPE()/ps,vertexBox2.t/ps,tpath2/ps,
				(siPMHit2.TimeFirstPE()-vertexBox2.t-tpath2)/ps))
		
		self.m.log(2,
			"dt =%7.2f ps, dtg= %7.2f ps, dpg =%7.2f ps, dtof= %7.2f ps,  "%(
				dt/ps,dtg/ps,dpg/ps,dtof/ps))
		

	############################################################		
	def Histos(self):
		"""
		book the histograms for the algo 
		"""
		self.TimeFirstPeBox1_histo_name = self.defineHisto("TimeFirstPEBox1",
										 "Time first PE Box1 (ps)",
										 1,[100],[-100],[400])

		self.TimeFirstPeBox2_histo_name = self.defineHisto("TimeFirstPEBox2",
										 "Time first PE Box2 (ps)",
										 1,[100],[-100],[400])

		self.TimeFirstPeBox1Box2_histo_name = self.defineHisto("TimeFirstPEBox1Box2",
										 "Time first PE Box1 vs Box2 (ps)",
										 2,[50,50],[-100,-100],[400,400])

		self.TimeScintLightBox1_histo_name = self.defineHisto("TimeScintLightBox1",
										 "Time scintillation ligth Box1 (ps)",
										 1,[50],[-100],[400])

		self.TimeScintLightBox2_histo_name = self.defineHisto("TimeScintLightBox2",
										 "Time scintillation ligth Box2 (ps)",
										 1,[50],[-100],[400])
		
		self.dpg_histo_name = self.defineHisto("DPG",
										 "Diff path scintillation light (ps)",
										 1,[50],[-500],[500])

		self.dtg_histo_name = self.defineHisto("DTG",
										 "Diff path gammas (ps)",
										 1,[50],[-500],[500])

		self.dt_histo_name = self.defineHisto("DT1stPE",
										 "Diff path 1st ps (ps)",
										 1,[50],[-500],[500])

		self.dtof_histo_name = self.defineHisto("DTOF",
										 " DTOF (ps)",
										 1,[50],[-500],[500])

		self.dtof2_histo_name = self.defineHisto("DTOF2",
										 " DTOF2 (ps)",
										 1,[80],[-200],[200])

		self.dtof3_histo_name = self.defineHisto("DTOF3",
										 " DTOF3 (ps)",
										 1,[80],[-100],[100])
		

	############################################################
	def defineHisto(self,histoName,histoTitle,histoType,nbin,xmin,xmax):

		histo_desc = histoName
		histo_name = self.alabel(histo_desc)
		if histoType == 1:
			self.hman.h1(histo_name, histo_desc,nbin[0],xmin[0],xmax[0])
		elif histoType == 2:
			self.hman.h2(histo_name, histo_desc,
				nbin[0],xmin[0],xmax[0],
				nbin[1],xmin[1],xmax[1])
		elif histoType == 3:
			self.hman.h3(histo_name, histo_desc,
				nbin[0],xmin[0],xmax[0],
				nbin[1],xmin[1],xmax[1],
				nbin[2],xmin[2],xmax[2])
		else:
			print "not implemented"
			sys.exit()
			
		self.hman.fetch(histo_name).GetXaxis().SetTitle(histoTitle)
		return histo_name
			
	
	