
def Trigger(self):
		"""
		Trigger the event 
		"""

		# #print "trigger"
		# ttr = -999
		
		# for ibin in range(NBINS):
		# 	pes = self.hsper.GetBinContent(ibin)
		# 	time = self.hsper.GetBinCenter(ibin)
		# 	#print "time (ps) pes =", time,pes
		# 	if pes > self.tLevel:
		# 		ttr = time
		# 		break

		# print "trigger time1 =",ttr
		# #wait()
		# hCTR1.Fill(ttr)

		# ttr = -999
		
		# if DEBUG > 1:
		# 	print self.STIME
		# 	print self.SPER

		# for ibin in xrange(len(self.STIME)):
		# 	time = self.STIME[ibin]
		# 	pes = self.SPER[ibin]
		# 	if DEBUG > 1:
		# 		print "time (ps) pes =", time,pes

		# 	if pes > self.tLevel:
		# 		ttr = time
		# 		break

		# if DEBUG > 0:
		# 	print "trigger time2 ps =",ttr
		# 	wait()
	
		# hCTR2.Fill(ttr)

# stime =[]

		# for tt in drange(0.,TAUMAX,TBINS):
		# 	stime.append(tt/ps)

		# self.STIME = np.array(stime)

#self.hsper = TH1F("hsper", "spe(au)", NBINS, 0.,(TAUMAX/ps))

		# if DEBUG > 0:
		# 	self.hSPER = TH1F("hSPER", "spe(au)", NBINS, 0.,(TAUMAX/ps))
		# 	self.c2 = TCanvas( 'c2', 'SPER', 200, 10, 600, 800 )

		# self.SPER = np.zeros(len(self.STIME))

# xx = -9999
		# yy = -9999
		# zz = -9999

		# if x0 == 0 or lx-x0 == 0:  #face x
		# 	return 0
		# if y0 == 0 or ly-y0 == 0:  #face y
		# 	return 1
		# if z0 == 0 or lz-z0 == 0:  #face z
		# 	return 2
		if -x0/tx >0 :
			D.append(-x0/tx)
			xx = 0
		elif (lx-x0)/tx >0 :
			D.append((lx-x0)/tx)
			xx = lx
		else:
			log.error(' negative distance to plane x: -x0/tx (mm) = %7.2f (lx-x0)/tx (mm) = %7.2f', 
				-(x0/tx)/mm,((lx-x0)/tx)/mm)
			sys.exit()

		if -y0/ty >0 :
			D.append(-y0/ty)
			yy = 0
		elif (ly-y0)/ty >=0 :
			D.append((ly-y0)/ty)
			yy = ly
		else:
			log.error(' negative distance to plane y: -y0/ty (mm) = %7.2f (ly-y0)/ty (mm) = %7.2f', 
				-(y0/ty)/mm,((ly-y0)/ty)/mm)
			sys.exit()

		if -z0/tz >0 :
			D.append(-z0/tz)
			zz = 0
			doi = photon.DOI()
		elif (lz-z0)/tz >0 :
			D.append((lz-z0)/tz)
			zz = lz
			doi = lz - photon.DOI()
		else:
			log.error(' negative distance to plane z: -z0/tz (mm) = %7.2f (lz-z0)/tz (mm) = %7.2f', 
				-(z0/tz)/mm,((lz-z0)/tz)/mm)
			sys.exit()

		if x == -9999: # generate random
			self.x = rnd.random()*self.WBox.X()
		elif x > 0 and x < self.WBox.X():
			self.x = x
		else:
			log.error('error: x outside box: x (cm) = %7.2f',x/cm)
			sys.exit()
		if y == -9999:
			self.y = rnd.random()*self.WBox.Y()
		elif y > 0 and y < self.WBox.Y():
			self.y = y
		else:
			log.error('error: x outside box: y (cm) = %7.2f',y/cm)
			sys.exit()
		if z == -9999:
			self.z = rnd.random()*self.WBox.Z()
		elif z > 0 and z < self.WBox.Z():
			self.z = z
		else:
			log.error('error: x outside box: z (cm) = %7.2f',z/cm)
			sys.exit()
		if tx == -9999:
			self.tx = rnd.uniform(-1.,1.)
		elif tx >= -1 and tx <= 1:
			self.tx = tx
		else:
			log.error('error: wrong tx = %7.2f',tx)
			sys.exit()
		if ty == -9999:
			self.ty = rnd.uniform(-1.,1.)
		elif ty >= -1 and ty <= 1:
			self.ty = ty
		else:
			log.error('error: wrong ty = %7.2f',ty)
			sys.exit()
		if tz == -9999:
			self.tz = rnd.uniform(-1.,1.)
		elif tz >= -1 and tz <= 1:
			self.tz = tz
		else:
			log.error('error: wrong tz = %7.2f',tz)
			sys.exit()

