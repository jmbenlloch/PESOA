from Centella.AAlgo import AAlgo
from Util import *
from Geometry import *
from TOF import *
from Scintillator import *
import random as rnd
import numpy


"""
This algorithm computes the Time map. 
"""


class DCR(AAlgo):

	############################################################
	def __init__(self, param=False, level=1, label="", **kargs):

		"""
		DCR Algorithm
		"""
		#self.m.log(1, 'Constructor()')

		### GENERAL STUFF
		self.name = 'DCR'
		
		AAlgo.__init__(self, param, level, self.name, 0, label, kargs)

    ### PARAMETERS
		
		self.debug = self.ints["Debug"]  #used to stop program at key break points

  		self.DCR = self.doubles["DCR"]  #Dark Count Rate in kHz

                self.waveforms = {} # Map {sensorID : waveform}, where waveform is a list [(bin,charge)]
 	
		if self.debug == 1:
			wait()


	############################################################		
	def initialize(self):

		self.m.log(1, 'Initialize()')
		
		### Defining histos
		self.Histos()

    ### Counters:
		self.numInputEvents = 0
		self.numOutputEvents = 0

	#Time Map instance
	#	self.sensorhits = TimeMap(numberOfBoxes = 2, dtmax=self.DTMAX)
		
		return


	############################################################
	def execute(self, event=""):

		self.m.log(2, 'Execute()')		
	
		self.numInputEvents += 1   

		#Add dark current to the waveform
                self.addDarkCurrent(event)
		
		if self.debug == 1:
			wait()	
		
		self.numOutputEvents += 1
		self.logman["USER"].gparam["waveforms"]=self.waveforms
		return True

############################################################
	def finalize(self):

#		self.m.log(1, 'Finalize()')

		self.m.log(1, 'Input  Events: ', self.numInputEvents)
		self.m.log(1, 'Output Events: ', self.numOutputEvents)

		return



############################################################
	def addDarkCurrent(self,event):
		"""
		Compute time maps
		"""

		sensorhitsTrue =  event.GetMCSensHits()
		self.m.log(3, " Add dark current: event has %d sensor hits = "%(
			len(sensorhitsTrue)))

		for hit in sensorhitsTrue:
			hid = hit.GetSensorID()
			waveform = hit.GetWaveform().GetData()
                        waveformDC = []

                        position = 0

                        lastBin = waveform[len(waveform)-1].first
                        print len(waveform)
                        print waveform[9].first

                        number = 0
                        for i in xrange(0,lastBin+1):
                            n = numpy.random.poisson(self.DCR*1000 * 5*10**(-12)) #Bins 5 ps
                            
                            print "i: %s, w: %s" % (i,waveform[position].first)

                            if i == waveform[position].first:
 #                               print "i: %s, w: %s" % (i,waveform[position].first)
                                number+=1
                                n = n + waveform[position].second
                                position+=1                            
                            if n > 0:
                                waveformDC.append((i,n))

                        self.waveforms[hid] = waveformDC

                        if len(waveformDC) != len(waveform):
                            print "!!! Noise added: %s - %s, number: %s" % (len(waveformDC),len(waveform),number)

                        for w in waveform:
                            print "Waveform: time bin: %s \tvalue: %s" % (w.first,w.second)

			self.m.log(5, " waveform for hit ID = %d: length = %d"%(
				hid,len(waveform)))
			
                        self.hman.fill(self.lastBinWaveform_histo_name,waveform[len(waveform)-1].first)
		

############################################################		
	def Histos(self):

            lastBinWaveform_histo_desc = "lastBinWaveform"
            self.lastBinWaveform_histo_name = self.alabel(lastBinWaveform_histo_desc)
            self.hman.h1(self.lastBinWaveform_histo_name, lastBinWaveform_histo_desc, 
                100, 0, 100000)
            self.hman.fetch(
            self.lastBinWaveform_histo_name).GetXaxis().SetTitle(
            "Last bin in the waveform")
