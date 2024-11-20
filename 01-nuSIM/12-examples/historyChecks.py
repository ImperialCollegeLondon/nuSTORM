'''

@file		historyChecks.py

	@brief		Provides functions for checking the eventHistory structure for nuSIM

	@author		Paul Kyberd

	Only plot muon neutrinos at the detector front face and not the whole front face plane
	@version	1.8
	@date		15 November 2024
	@author		Paul Kyberd		

	Add cuts for x, y, xp
	@version	1.7
	@date		14 November 2024
	@author		Paul Kyberd		

	Add reading in a cuts dictioanary using json
	@version	1.6
	@date		13 November 2024
	@author		Paul Kyberd		

	Add a muon decay station
	@version	1.5
	@date		13 November 2024
	@author		Paul Kyberd		

	Add a numDetector station
	@version	1.4
	@date		27 October 2024
	@author		Paul Kyberd	

	Add a piFlashNu station
	@version	1.3
	@date		26 October 2024
	@author		Paul Kyberd	

	Add a pdgcode plot
	@version	1.2
	@date		25 October 2024
	@author		Paul Kyberd

	Add a parameter in instantiation to set cut on weight or not - doesn't work - do by hand
	@version	1.1
	@date		25 October 2024
	@author		Paul Kyberd

	To instantiate pass the pointer to the histomanager object
	@version	1.0
	@date		24 October 2024
	@author		Paul Kyberd

'''
import histoManager as histoManager
import json

class historyChecks:


	def __init__(self, hm):
		print("instantiating")
		self._hm = hm
		self._targetPlots = False
		self._prodStrghtStartPlots = False
		self._prodStrghtEndPlots = False
		self._numuDetPlots = False
		self._piFlashNuPlots = False
		self._muonDecayPlots = False
		self._plotSet = False
		self._particleSelect = False


	def __str__(self):
		instanceStr = "string from historyChecks" 
		return instanceStr 

	def __repr__(self):
		instanceStr = "string from historyChecks" 
		return instanceStr 

	def createTargetPlots(self, wtCutFlg, cutsDict, particleID):
		self._wtCutFlg = wtCutFlg
		self._particleID = particleID
		if self._plotSet:
			exit("Trying to create target plots - other station plots already created")
		self._targetPlots = True
		station = "Target"
		self.createPlots(station, cutsDict)

	def createPStStartPlots(self, wtCutFlg, cutsDict, particleID):
		self._wtCutFlg = wtCutFlg
		self._particleID = particleID
		if self._plotSet:
			exit("Trying to create production straight start plots - other station plots already created")
		self._prodStrghtStartPlots = True
		station = "productionStraight "
		self.createPlots(station, cutsDict)

	def createPSEndPlots(self, wtCutFlg, cutsDict, particleID):
		self._wtCutFlg = wtCutFlg
		self._particleID = particleID
		if self._plotSet:
			exit("Trying to create production straight end plots - other station plots already created")
		self._prodStrghtEndPlots = True
		station = "prodStraightEnd"
		if self._particleID == 0: 
			print ("no particle selected ")
			self._particleSelectFlg = False
		else:
			self._particleID = particleID
			self._particleSelectFlg = True

		self.createPlots(station, cutsDict)

	def createnumuDetectorPlots(self, wtCutFlg, cutsDict, particleID):
		self._wtCutFlg = wtCutFlg
		self._particleID = particleID
		if self._plotSet:
			exit("Trying to create numu detector plots - other station plots already created")
		self._numuDetPlots = True
		station = "numuDetector"
		if self._particleID == 0: 
			print ("no particle selected ")
			self._particleSelectFlg = False
		else:
			self._particleID = particleID
			self._particleSelect = True
		self.createPlots(station, cutsDict)

	def createPiFlashNuPlots(self, wtCutFlg, cutsDict, particleID):
		self._wtCutFlg = wtCutFlg
		self._particleID = particleID
		if self._plotSet:
			exit("Trying to create pi flash nu plots - other station plots already created")
		self._piFlashNuPlots = True
		station = "piFlashNu"
		if self._particleID == 0: 
			print ("no particle selected ")
			self._particleSelectFlg = False
		else:
			self._particleID = particleID
			self._particleSelect = True
		self.createPlots(station, cutsDict)

	def createmuonDecayPlots(self, wtCutFlg, cutsDict, particleID):
		self._wtCutFlg = wtCutFlg
		self._particleID = particleID
		if self._plotSet:
			exit("Trying to create muon decay plots - other station plots already created")
		self._muonDecayPlots = True
		station = "muonDecay"
		if self._particleID == 0: 
			print ("no particle selected ")
			self._particleSelectFlg = False
		else:
			self._particleID = particleID
			self._particleSelect = True
		self.createPlots(station, cutsDict)


	def createPlots(self, station, cutsDict):
		print ("in creating plots for station ", station)

		#	Read the cuts file
		with open(cutsDict) as pD:
			self._cutsInfo = json.load(pD)


		hBins  = 100
		hLow = 0.0 
		hUp  = 0.0

#	set any lower or upper values found in the file
#	x
		try:
			hLowX = self._cutsInfo['xLower'][station]
		except KeyError:
			hLowX = 0.0
		try:
			hUpX = self._cutsInfo['xHigher'][station]
		except KeyError:
			hUpX = 0.0
#	y
		try:
			hLowY = self._cutsInfo['yLower'][station]
		except KeyError:
			hLowY = 0.0
		try:
			hUpY = self._cutsInfo['yHigher'][station]
		except KeyError:
			hUpY = 0.0
#	xp
		try:
			hLowXP = self._cutsInfo['xpLower'][station]
		except KeyError:
			hLowXP = 0.0
		try:
			hUpXP = self._cutsInfo['xpHigher'][station]
		except KeyError:
			hUpXP = 0.0
#	yp
		try:
			hLowYP = self._cutsInfo['ypLower'][station]
		except KeyError:
			hLowYP = 0.0
		try:
			hUpYP = self._cutsInfo['ypHigher'][station]
		except KeyError:
			hUpYP = 0.0
#	px
		try:
			hLowPx = self._cutsInfo['pxLower'][station]
		except KeyError:
			hLowPx = 0.0
		try:
			hUpPx = self._cutsInfo['pxHigher'][station]
		except KeyError:
			hUpPx = 0.0


		print ("hLowXP is ", hLowXP, "    hUpXP is ", hUpXP)

		hTitle = station + ": x"
		self.htx = self._hm.book(hTitle, hBins, hLowX, hUpX)
		hTitle = station + ": y"
		self.hty = self._hm.book(hTitle, hBins, hLowY, hUpY)
		hTitle = station + ": z"
		self.htz = self._hm.book(hTitle, hBins, hLow, hUp)
		hTitle = station + ": s"
		self.hts = self._hm.book(hTitle, hBins, hLow, hUp)
		hTitle = station + ": t"
		self.htt = self._hm.book(hTitle, hBins, hLow, hUp)

		hTitle = station + ": Px"
		self.htpx = self._hm.book(hTitle, hBins, hLowPx, hUpPx)
		hTitle = station + ": Py"
		self.htpy = self._hm.book(hTitle, hBins, hLow, hUp)
		hTitle = station + ": Pz"
		self.htpz = self._hm.book(hTitle, hBins, hLow, hUp)

		hTitle = station + ": Xp"
		self.htxp = self._hm.book(hTitle, hBins, hLowXP, hUpXP)
		hTitle = station + ": Yp"
		self.htyp = self._hm.book(hTitle, hBins, hLow, hUp)

		hTitle = station + ": p"
		self.htp = self._hm.book(hTitle, hBins, hLow, hUp)

		hTitle = station + ": x v y"
		self.htxy = self._hm.book2(hTitle, hBins, hLowX, hUpX, hBins, hLowY, hUpY)
		hTitle = station + ": Xp v Yp"
		self.htxpyp = self._hm.book2(hTitle, hBins, hLowXP, hUpXP, hBins, hLowYP, hUpYP)
		hTitle = station + ": x v Xp"
		self.htxxp = self._hm.book2(hTitle, hBins, hLowX, hUpX, hBins, hLowXP, hUpXP)
		hTitle = station + ": y v Yp"
		self.htyyp = self._hm.book2(hTitle, hBins, hLowY, hUpY, hBins, hLow, hUp)

		hTitle = station + ": weight"
		self.hwt = self._hm.book(hTitle, hBins, hLow, hUp)

		hTitle = station + ": PDG code"
		self.hPDG = self._hm.book(hTitle, hBins, hLow, hUp)


		self._plots = True

	def fillTargetPlots(self, particle):
		if not self._targetPlots:
			exit("Trying to fill a plot for the wrong station")
		self.fillPlots(particle)

	def fillPSStartPlots(self, particle):
		if not self._prodStrghtStartPlots:
			exit("Trying to fill a plot for the wrong station")
		self.fillPlots(particle)

	def fillPSEndPlots(self, particle):
		if not self._prodStrghtEndPlots:
			exit("Trying to fill a plot for the wrong station")
		self.fillPlots(particle)

	def fillnumuDetectorPlots(self, particle):
		if not self._numuDetPlots:
			exit("Trying to fill a plot numuDetector, the wrong station")
		#	Onlt plot hits on the detector - not flux through the detector plane
		if particle.weight() > 1:
			self.fillPlots(particle)

	def fillPiFlashNuPlots(self, particle):
		if not self._piFlashNuPlots:
			exit("Trying to fill a plot pi Flash Nu, the wrong station")
		self.fillPlots(particle)

	def fillmuonDecayPlots(self, particle):
		if not self._muonDecayPlots:
			exit("Trying to fill a plot muon Decay, the wrong station")
		self.fillPlots(particle)


	def fillPlots(self, particle):

#	Only try to fill histograms if they have been initialised	
		if not self._plots:
			print ("plots not initialised")
			return

#
		weight = particle.weight()
		PDG = particle.pdgCode()

		fillFlag = False

		if not self._particleSelect:
			fillFlag = True
		elif self._particleID == PDG:
			fillFlag = True

#		if self._wtCutFlg == True and weight > 1: 
#			fillFlag = True
#		if not self._wtCutFlg: 
#			fillFlag = True

#		fillFlag = True
		if weight > 1 and fillFlag:
			fillFlag = True
		else:
			fillFlag = False

		if fillFlag:
			self.htx.Fill(particle.x())
			self.hty.Fill(particle.y())
			self.htz.Fill(particle.z())
			self.hts.Fill(particle.s())
			self.htt.Fill(particle.t())
			self.htpx.Fill(particle.px())
			self.htpy.Fill(particle.py())
			self.htpz.Fill(particle.pz())
			self.htxp.Fill(particle.xp())
			self.htyp.Fill(particle.yp())
			self.htp.Fill(particle.p()[0])

			self.htxy.Fill(particle.x(), particle.y())
			self.htxxp.Fill(particle.x(), particle.xp())
			self.htyyp.Fill(particle.y(), particle.yp())
			self.htxpyp.Fill(particle.xp(), particle.yp())

			self.hwt.Fill(weight)
			self.hPDG.Fill(particle.pdgCode())
		
		return


if __name__ == "__main__" :
	inst = historyChecks()
	print(inst)
