#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Class: eventHistory:
====================

    @file eventHistory.py

    @brief   creates a history of a nuStorm event

    @author  Paul Kyberd

    @version     1.1
    @date        11 November 2021

Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__: Initialisation of
      __repr__: One liner with call.
      __str__ : Summary

  Get/set methods:

   	makeHistory()		: Create an empty structure (all pi+) to be filled by the user with add particle
		getEntries()		: total number of entries in the root file (for reading)
		readNext()			: read next entry in the root file
		addParticle(location, par): add particle par at location to the history
		findParticle(location, par): return particle par at location from the history

  Generation, calcution methods and utilities:

   	outFile(fileName)	: Name and open a file for output
   	inFile(fileName)	: Name and open a file for input
		readNext()			: read next entry in the root file
		addParticle(location, par): add particle par at location to the history
		findParticle(location, par): return particle par at location from the history




Version 1.2 									29/04/2022
Add a point at the end of the production straight

Version 1.1										11/11/2021
After writing out

Version 1.0										16/09/2021
Initial version

'''

# Generic Python imports
import sys, os

import array
from copy import deepcopy
# numpy
import numpy as np

# pyRoot imports
import ROOT
from ROOT import TFile, TNtuple, gROOT, TROOT, TTree, AddressOf, addressof

# nuStorm imports
import particle as particle
import MuonConst as MuonConst
import PionConst as PionConst
import NeutrinoEventInstance as nuEvtInst

gROOT.ProcessLine(
"struct target {\
   Int_t				runNumber;\
   Int_t				eventNumber;\
   Int_t				pdgCode;\
   Float_t			x;\
   Float_t			y;\
   Float_t			z;\
   Float_t			s;\
   Float_t			px;\
   Float_t			py;\
   Float_t			pz;\
   Float_t			t;\
   Float_t			eventWeight;\
   Float_t			mass;\
};" );

gROOT.ProcessLine(
"struct productionStraight {\
   Int_t				runNumber;\
   Int_t				eventNumber;\
   Int_t				pdgCode;\
   Float_t			x;\
   Float_t			y;\
   Float_t			z;\
   Float_t			s;\
   Float_t			px;\
   Float_t			py;\
   Float_t			pz;\
   Float_t			t;\
   Float_t			eventWeight;\
   Float_t			mass;\
};" );

gROOT.ProcessLine(
"struct prodStraightEnd {\
   Int_t				runNumber;\
   Int_t				eventNumber;\
   Int_t				pdgCode;\
   Float_t			x;\
   Float_t			y;\
   Float_t			z;\
   Float_t			s;\
   Float_t			px;\
   Float_t			py;\
   Float_t			pz;\
   Float_t			t;\
   Float_t			eventWeight;\
   Float_t			mass;\
};" );


gROOT.ProcessLine(
"struct pionDecay {\
   Int_t				runNumber;\
   Int_t				eventNumber;\
   Int_t				pdgCode;\
   Float_t			x;\
   Float_t			y;\
   Float_t			z;\
   Float_t			s;\
   Float_t			px;\
   Float_t			py;\
   Float_t			pz;\
   Float_t			t;\
   Float_t			eventWeight;\
   Float_t			mass;\
};" );

gROOT.ProcessLine(
"struct muonProduction {\
   Int_t				runNumber;\
   Int_t				eventNumber;\
   Int_t				pdgCode;\
   Float_t			x;\
   Float_t			y;\
   Float_t			z;\
   Float_t			s;\
   Float_t			px;\
   Float_t			py;\
   Float_t			pz;\
   Float_t			t;\
   Float_t			eventWeight;\
   Float_t			mass;\
};" );

gROOT.ProcessLine(
"struct piFlashNu {\
   Int_t				runNumber;\
   Int_t				eventNumber;\
   Int_t				pdgCode;\
   Float_t			x;\
   Float_t			y;\
   Float_t			z;\
   Float_t			s;\
   Float_t			px;\
   Float_t			py;\
   Float_t			pz;\
   Float_t			t;\
   Float_t			eventWeight;\
   Float_t			mass;\
};" );

gROOT.ProcessLine(
"struct muonDecay {\
   Int_t				runNumber;\
   Int_t				eventNumber;\
   Int_t				pdgCode;\
   Float_t			x;\
   Float_t			y;\
   Float_t			z;\
   Float_t			s;\
   Float_t			px;\
   Float_t			py;\
   Float_t			pz;\
   Float_t			t;\
   Float_t			eventWeight;\
   Float_t			mass;\
};" );

gROOT.ProcessLine(
"struct eProduction {\
   Int_t				runNumber;\
   Int_t				eventNumber;\
   Int_t				pdgCode;\
   Float_t			x;\
   Float_t			y;\
   Float_t			z;\
   Float_t			s;\
   Float_t			px;\
   Float_t			py;\
   Float_t			pz;\
   Float_t			t;\
   Float_t			eventWeight;\
   Float_t			mass;\
};" );

gROOT.ProcessLine(
"struct numuProduction {\
   Int_t				runNumber;\
   Int_t				eventNumber;\
   Int_t				pdgCode;\
   Float_t			x;\
   Float_t			y;\
   Float_t			z;\
   Float_t			s;\
   Float_t			px;\
   Float_t			py;\
   Float_t			pz;\
   Float_t			t;\
   Float_t			eventWeight;\
   Float_t			mass;\
};" );

gROOT.ProcessLine(
"struct nueProduction {\
   Int_t				runNumber;\
   Int_t				eventNumber;\
   Int_t				pdgCode;\
   Float_t			x;\
   Float_t			y;\
   Float_t			z;\
   Float_t			s;\
   Float_t			px;\
   Float_t			py;\
   Float_t			pz;\
   Float_t			t;\
   Float_t			eventWeight;\
   Float_t			mass;\
};" );

gROOT.ProcessLine(
"struct numuDetector {\
   Int_t				runNumber;\
   Int_t				eventNumber;\
   Int_t				pdgCode;\
   Float_t			x;\
   Float_t			y;\
   Float_t			z;\
   Float_t			s;\
   Float_t			px;\
   Float_t			py;\
   Float_t			pz;\
   Float_t			t;\
   Float_t			eventWeight;\
   Float_t			mass;\
};" );

gROOT.ProcessLine(
"struct nueDetector {\
   Int_t				runNumber;\
   Int_t				eventNumber;\
   Int_t				pdgCode;\
   Float_t			x;\
   Float_t			y;\
   Float_t			z;\
   Float_t			s;\
   Float_t			px;\
   Float_t			py;\
   Float_t			pz;\
   Float_t			t;\
   Float_t			eventWeight;\
   Float_t			mass;\
};" );


class eventHistory:
	Version = 1.0
	__Validated__ = False

# built in methods
	def __init__(self ):
		self._outputFilename = "null"
		self._inputFilename = "null"
		self._eHTree = "null"
		self._runNum = -1
		self._particles = [None]*12
		self.evTree = TTree('eventHistory', 'nuStorm Event')
		self._entryPnt = 0

	def __repr__(self):
		return "a history of a nuSTORM event"

	def __str__(self):
		return "eventHistory: run Number %i, outputFilename = %s" \
		         % \
               (self._runNum, self._outputFilename)


#	methods
#  Create an empty structure to be filled by the user with add particle - at present the
#  the particles in the empty structure are all pi+
	def makeHistory(self):
#  make sure the data structure is complete
		testParticle = particle.particle(-1, -1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, "pi+")
		self.addParticle("target", testParticle)
		self.addParticle("productionStraight", testParticle)
		self.addParticle("prodStraightEnd", testParticle)
		self.addParticle("pionDecay", testParticle)
		self.addParticle("muonProduction", testParticle)
		self.addParticle("piFlashNu", testParticle)
		self.addParticle("muonDecay", testParticle)
		self.addParticle("eProduction", testParticle)
		self.addParticle("numuProduction", testParticle)
		self.addParticle("nueProduction", testParticle)
		self.addParticle("numuDetector", testParticle)
		self.addParticle("nueDetector", testParticle)

#  Name and open a file for output
	def outFile(self, fileName ):
		self._outputFilename = fileName
		self._outTFile = TFile( self._outputFilename, 'RECREATE', 'event History' )

#  Name and open a file for input
	def inFile(self, fileName ):
		self._inputFilename = fileName
		self._inTFile = TFile( self._inputFilename, 'READ', 'event History' )

		self._eHTree = self._inTFile.Get("eventHistory")
		self.target = ROOT.target()
		self._eHTree.SetBranchAddress('target', self.target)
		self.productionStraight = ROOT.productionStraight()
		self._eHTree.SetBranchAddress('productionStraight', self.productionStraight)
		self.prodStraightEnd = ROOT.prodStraightEnd()
		self._eHTree.SetBranchAddress('prodStraightEnd', self.prodStraightEnd)
		self.pionDecay = ROOT.pionDecay()
		self._eHTree.SetBranchAddress("pionDecay", self.pionDecay)
		self.muonProduction = ROOT.muonProduction()
		self._eHTree.SetBranchAddress("muonProduction", self.muonProduction)
		self.piFlashProd = ROOT.piFlashNu()
		self._eHTree.SetBranchAddress("piFlashNu", self.piFlashProd)
		self.muonDecay = ROOT.muonDecay()
		self._eHTree.SetBranchAddress("muonDecay", self.muonDecay)
		self.eProduction = ROOT.eProduction()
		self._eHTree.SetBranchAddress("eProduction", self.eProduction)
		self.numuProduction = ROOT.numuProduction()
		self._eHTree.SetBranchAddress("numuProduction", self.numuProduction)
		self.nueProduction = ROOT.nueProduction()
		self._eHTree.SetBranchAddress("nueProduction", self.nueProduction)
		self.numuDetector = ROOT.numuDetector()
		self._eHTree.SetBranchAddress("numuDetector", self.numuDetector)
		self.nueDetector = ROOT.nueDetector()
		self._eHTree.SetBranchAddress("nueDetector", self.nueDetector)

	def getEntries(self):
		return self._eHTree.GetEntries()

	def readNext(self):
#		eHTree.Show(0)
#		print ("number of entries is ", eHTree.GetEntries())
#		eHTree.Show(1)

#		for b in eHTree.GetListOfBranches():
#			print("branch", b.GetName())
#			print(b)

#		for pnt in range(2):
#			self._eHTree.GetEntry(pnt)
#			print ("at target x ", self.target.x)
#			print ("production Straight x ", self.productionStraight.x)
#			print ("pion Decay x ", self.pionDecay.x)
#			print (self.pionDecay.runNumber)
#			print (self.pionDecay.eventNumber)

		self._eHTree.GetEntry(self._entryPnt)

		pTar = particle.particle(self.target.runNumber, self.target.eventNumber, self.target.s,
			self.target.x, self.target.y, self.target.z, self.target.px, self.target.py, self.target.pz,
			self.target.t, self.target.eventWeight, self.target.pdgCode)
		pPs = particle.particle(self.productionStraight.runNumber, self.productionStraight.eventNumber, self.productionStraight.s,
			self.productionStraight.x, self.productionStraight.y, self.productionStraight.z, self.productionStraight.px, self.productionStraight.py,
			self.productionStraight.pz, self.productionStraight.t, self.productionStraight.eventWeight, self.productionStraight.pdgCode)
		pPsEnd = particle.particle(self.prodStraightEnd.runNumber, self.prodStraightEnd.eventNumber, self.prodStraightEnd.s,
			self.prodStraightEnd.x, self.prodStraightEnd.y, self.prodStraightEnd.z, self.prodStraightEnd.px, self.prodStraightEnd.py,
			self.prodStraightEnd.pz, self.prodStraightEnd.t, self.prodStraightEnd.eventWeight, self.prodStraightEnd.pdgCode)
		pPiDcy = particle.particle(self.pionDecay.runNumber, self.pionDecay.eventNumber, self.pionDecay.s,
			self.pionDecay.x, self.pionDecay.y, self.pionDecay.z, self.pionDecay.px, self.pionDecay.py, self.pionDecay.pz,
			self.pionDecay.t, self.pionDecay.eventWeight, self.pionDecay.pdgCode)
		pMuPrd = particle.particle(self.muonProduction.runNumber, self.muonProduction.eventNumber, self.muonProduction.s,
			self.muonProduction.x, self.muonProduction.y, self.muonProduction.z, self.muonProduction.px, self.muonProduction.py, self.muonProduction.pz,
			self.muonProduction.t, self.muonProduction.eventWeight, self.muonProduction.pdgCode)
		pFlsNu = particle.particle(self.piFlashProd.runNumber, self.piFlashProd.eventNumber, self.piFlashProd.s,
			self.piFlashProd.x, self.piFlashProd.y, self.piFlashProd.z, self.piFlashProd.px, self.piFlashProd.py, self.piFlashProd.pz,
			self.piFlashProd.t, self.piFlashProd.eventWeight, self.piFlashProd.pdgCode)
		pMuDcy = particle.particle(self.muonDecay.runNumber, self.muonDecay.eventNumber, self.muonDecay.s,
			self.muonDecay.x, self.muonDecay.y, self.muonDecay.z, self.muonDecay.px, self.muonDecay.py, self.muonDecay.pz,
			self.muonDecay.t, self.muonDecay.eventWeight, self.muonDecay.pdgCode)
		pEPrd = particle.particle(self.eProduction.runNumber, self.eProduction.eventNumber, self.eProduction.s,
			self.eProduction.x, self.eProduction.y, self.eProduction.z, self.eProduction.px, self.eProduction.py, self.eProduction.pz,
			self.eProduction.t, self.eProduction.eventWeight, self.eProduction.pdgCode)
		pNumuPrd = particle.particle(self.numuProduction.runNumber, self.numuProduction.eventNumber, self.numuProduction.s,
			self.numuProduction.x, self.numuProduction.y, self.numuProduction.z, self.numuProduction.px, self.numuProduction.py, self.numuProduction.pz,
			self.numuProduction.t, self.numuProduction.eventWeight, self.numuProduction.pdgCode)
		pNuePrd = particle.particle(self.numuProduction.runNumber, self.numuProduction.eventNumber, self.nueProduction.s,
			self.nueProduction.x, self.nueProduction.y, self.nueProduction.z, self.nueProduction.px, self.nueProduction.py, self.nueProduction.pz,
			self.nueProduction.t, self.nueProduction.eventWeight, self.nueProduction.pdgCode)
		pNumuDet = particle.particle(self.numuDetector.runNumber, self.numuDetector.eventNumber, self.numuDetector.s,
			self.numuDetector.x, self.numuDetector.y, self.numuDetector.z, self.numuDetector.px, self.numuDetector.py, self.numuDetector.pz,
			self.numuDetector.t, self.numuDetector.eventWeight, self.numuDetector.pdgCode)
		pNueDet = particle.particle(self.nueDetector.runNumber, self.nueDetector.eventNumber, self.nueDetector.s,
			self.nueDetector.x, self.nueDetector.y, self.nueDetector.z, self.nueDetector.px, self.nueDetector.py, self.nueDetector.pz,
			self.nueDetector.t, self.nueDetector.eventWeight, self.nueDetector.pdgCode)

		self.addParticle('target', pTar)
		self.addParticle('productionStraight', pPs)
		self.addParticle('prodStraightEnd', pPsEnd)
		self.addParticle('pionDecay', pPiDcy)
		self.addParticle('muonProduction', pMuPrd)
		self.addParticle('piFlashNu', pFlsNu)
		self.addParticle('muonDecay', pMuDcy)
		self.addParticle('eProduction', pEPrd)
		self.addParticle('numuProduction', pNumuPrd)
		self.addParticle('nueProduction', pNuePrd)
		self.addParticle('numuDetector', pNumuDet)
		self.addParticle('nueDetector', pNueDet)
		self._entryPnt = self._entryPnt + 1

	def rootStructure(self):
# Define the root data structure
		self.target = ROOT.target()
		self.evTree.Branch('target', self.target, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
		self.productionStraight = ROOT.productionStraight()
		self.evTree.Branch('productionStraight', self.productionStraight, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
		self.prodStraightEnd = ROOT.prodStraightEnd()
		self.evTree.Branch('prodStraightEnd', self.prodStraightEnd, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
		self.pionDecay = ROOT.pionDecay()
		self.evTree.Branch('pionDecay', self.pionDecay, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
		self.muonProduction = ROOT.muonProduction()
		self.evTree.Branch('muonProduction', self.muonProduction, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
		self.piFlashNu = ROOT.piFlashNu()
		self.evTree.Branch('piFlashNu', self.piFlashNu, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
		self.muonDecay = ROOT.muonDecay()
		self.evTree.Branch('muonDecay', self.muonDecay, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
		self.eProduction = ROOT.eProduction()
		self.evTree.Branch('eProduction', self.eProduction, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
		self.numuProduction = ROOT.numuProduction()
		self.evTree.Branch('numuProduction', self.numuProduction, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
		self.nueProduction = ROOT.nueProduction()
		self.evTree.Branch('nueProduction', self.nueProduction, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
		self.numuDetector = ROOT.numuDetector()
		self.evTree.Branch('numuDetector', self.numuDetector, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
		self.nueDetector = ROOT.nueDetector()
		self.evTree.Branch('nueDetector', self.nueDetector, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
# and create the history array
		self.makeHistory()

#  Close the output file
	def outFileClose(self ):
		self._outTFile.Close()

	def fill(self):

# target
		partPnt = 0
		self.target.runNumber = self._particles[0].run()
		self.target.eventNumber = self._particles[0].event()
		self.target.pdgCode = self._particles[0].pdgCode()
		self.target.x = self._particles[0].x()
		self.target.y = self._particles[0].y()
		self.target.z = self._particles[0].z()
		self.target.s = self._particles[0].s()
		self.target.px = self._particles[0].p()[1][0]
		self.target.py = self._particles[0].p()[1][1]
		self.target.pz = self._particles[0].p()[1][2]
		self.target.t = self._particles[0].t()
		self.target.eventWeight = self._particles[0].weight()
		self.target.mass = self._particles[0].mass()

# production straight (start)
		partPnt = partPnt + 1
		self.productionStraight.runNumber = self._particles[partPnt].run()
		self.productionStraight.eventNumber = self._particles[partPnt].event()
		self.productionStraight.pdgCode = self._particles[partPnt].pdgCode()
		self.productionStraight.x = self._particles[partPnt].x()
		self.productionStraight.y = self._particles[partPnt].y()
		self.productionStraight.z = self._particles[partPnt].z()
		self.productionStraight.s = self._particles[partPnt].s()
		self.productionStraight.px = self._particles[partPnt].p()[1][0]
		self.productionStraight.py = self._particles[partPnt].p()[1][1]
		self.productionStraight.pz = self._particles[partPnt].p()[1][2]
		self.productionStraight.t = self._particles[partPnt].t()
		self.productionStraight.eventWeight = self._particles[partPnt].weight()
		self.productionStraight.mass = self._particles[partPnt].mass()

# production straight end
		partPnt = partPnt + 1
		self.prodStraightEnd.runNumber = self._particles[partPnt].run()
		self.prodStraightEnd.eventNumber = self._particles[partPnt].event()
		self.prodStraightEnd.pdgCode = self._particles[partPnt].pdgCode()
		self.prodStraightEnd.x = self._particles[partPnt].x()
		self.prodStraightEnd.y = self._particles[partPnt].y()
		self.prodStraightEnd.z = self._particles[partPnt].z()
		self.prodStraightEnd.s = self._particles[partPnt].s()
		self.prodStraightEnd.px = self._particles[partPnt].p()[1][0]
		self.prodStraightEnd.py = self._particles[partPnt].p()[1][1]
		self.prodStraightEnd.pz = self._particles[partPnt].p()[1][2]
		self.prodStraightEnd.t = self._particles[partPnt].t()
		self.prodStraightEnd.eventWeight = self._particles[partPnt].weight()
		self.prodStraightEnd.mass = self._particles[partPnt].mass()

# pion decay
		partPnt = partPnt + 1
		self.pionDecay.runNumber = self._particles[partPnt].run()
		self.pionDecay.eventNumber = self._particles[partPnt].event()
		self.pionDecay.pdgCode = self._particles[partPnt].pdgCode()
		self.pionDecay.x = self._particles[partPnt].x()
		self.pionDecay.y = self._particles[partPnt].y()
		self.pionDecay.z = self._particles[partPnt].z()
		self.pionDecay.s = self._particles[partPnt].s()
		self.pionDecay.px = self._particles[partPnt].p()[1][0]
		self.pionDecay.py = self._particles[partPnt].p()[1][1]
		self.pionDecay.pz = self._particles[partPnt].p()[1][2]
		self.pionDecay.t = self._particles[partPnt].t()
		self.pionDecay.eventWeight = self._particles[partPnt].weight()
		self.pionDecay.mass = self._particles[partPnt].mass()

# muon production
		partPnt = partPnt + 1
		self.muonProduction.runNumber = self._particles[partPnt].run()
		self.muonProduction.eventNumber = self._particles[partPnt].event()
		self.muonProduction.pdgCode = self._particles[partPnt].pdgCode()
		self.muonProduction.x = self._particles[partPnt].x()
		self.muonProduction.y = self._particles[partPnt].y()
		self.muonProduction.z = self._particles[partPnt].z()
		self.muonProduction.s = self._particles[partPnt].s()
		self.muonProduction.px = self._particles[partPnt].p()[1][0]
		self.muonProduction.py = self._particles[partPnt].p()[1][1]
		self.muonProduction.pz = self._particles[partPnt].p()[1][2]
		self.muonProduction.t = self._particles[partPnt].t()
		self.muonProduction.eventWeight = self._particles[partPnt].weight()
		self.muonProduction.mass = self._particles[partPnt].mass()

# neutrino flash production
		partPnt = partPnt + 1
		self.piFlashNu.runNumber = self._particles[partPnt].run()
		self.piFlashNu.eventNumber = self._particles[partPnt].event()
		self.piFlashNu.pdgCode = self._particles[partPnt].pdgCode()
		self.piFlashNu.x = self._particles[partPnt].x()
		self.piFlashNu.y = self._particles[partPnt].y()
		self.piFlashNu.z = self._particles[partPnt].z()
		self.piFlashNu.s = self._particles[partPnt].s()
		self.piFlashNu.px = self._particles[partPnt].p()[1][0]
		self.piFlashNu.py = self._particles[partPnt].p()[1][1]
		self.piFlashNu.pz = self._particles[partPnt].p()[1][2]
		self.piFlashNu.t = self._particles[partPnt].t()
		self.piFlashNu.eventWeight = self._particles[partPnt].weight()
		self.piFlashNu.mass = self._particles[partPnt].mass()

# muon decay
		partPnt = partPnt + 1
		self.muonDecay.runNumber = self._particles[partPnt].run()
		self.muonDecay.eventNumber = self._particles[partPnt].event()
		self.muonDecay.pdgCode = self._particles[partPnt].pdgCode()
		self.muonDecay.x = self._particles[partPnt].x()
		self.muonDecay.y = self._particles[partPnt].y()
		self.muonDecay.z = self._particles[partPnt].z()
		self.muonDecay.s = self._particles[partPnt].s()
		self.muonDecay.px = self._particles[partPnt].p()[1][0]
		self.muonDecay.py = self._particles[partPnt].p()[1][1]
		self.muonDecay.pz = self._particles[partPnt].p()[1][2]
		self.muonDecay.t = self._particles[partPnt].t()
		self.muonDecay.eventWeight = self._particles[partPnt].weight()
		self.muonDecay.mass = self._particles[partPnt].mass()

# e production
		partPnt = partPnt + 1
		self.eProduction.runNumber = self._particles[partPnt].run()
		self.eProduction.eventNumber = self._particles[partPnt].event()
		self.eProduction.pdgCode = self._particles[partPnt].pdgCode()
		self.eProduction.x = self._particles[partPnt].x()
		self.eProduction.y = self._particles[partPnt].y()
		self.eProduction.z = self._particles[partPnt].z()
		self.eProduction.s = self._particles[partPnt].s()
		self.eProduction.px = self._particles[partPnt].p()[1][0]
		self.eProduction.py = self._particles[partPnt].p()[1][1]
		self.eProduction.pz = self._particles[partPnt].p()[1][2]
		self.eProduction.t = self._particles[partPnt].t()
		self.eProduction.eventWeight = self._particles[partPnt].weight()
		self.eProduction.mass = self._particles[partPnt].mass()

# numu production
		partPnt = partPnt + 1
		self.numuProduction.runNumber = self._particles[partPnt].run()
		self.numuProduction.eventNumber = self._particles[partPnt].event()
		self.numuProduction.pdgCode = self._particles[partPnt].pdgCode()
		self.numuProduction.x = self._particles[partPnt].x()
		self.numuProduction.y = self._particles[partPnt].y()
		self.numuProduction.z = self._particles[partPnt].z()
		self.numuProduction.s = self._particles[partPnt].s()
		self.numuProduction.px = self._particles[partPnt].p()[1][0]
		self.numuProduction.py = self._particles[partPnt].p()[1][1]
		self.numuProduction.pz = self._particles[partPnt].p()[1][2]
		self.numuProduction.t = self._particles[partPnt].t()
		self.numuProduction.eventWeight = self._particles[partPnt].weight()
		self.numuProduction.mass = self._particles[partPnt].mass()

# nue production
		partPnt = partPnt + 1
		self.nueProduction.runNumber = self._particles[partPnt].run()
		self.nueProduction.eventNumber = self._particles[partPnt].event()
		self.nueProduction.pdgCode = self._particles[partPnt].pdgCode()
		self.nueProduction.x = self._particles[partPnt].x()
		self.nueProduction.y = self._particles[partPnt].y()
		self.nueProduction.z = self._particles[partPnt].z()
		self.nueProduction.s = self._particles[partPnt].s()
		self.nueProduction.px = self._particles[partPnt].p()[1][0]
		self.nueProduction.py = self._particles[partPnt].p()[1][1]
		self.nueProduction.pz = self._particles[partPnt].p()[1][2]
		self.nueProduction.t = self._particles[partPnt].t()
		self.nueProduction.eventWeight = self._particles[partPnt].weight()
		self.nueProduction.mass = self._particles[partPnt].mass()

# numu Detector
		partPnt = partPnt + 1
		self.numuDetector.runNumber = self._particles[partPnt].run()
		self.numuDetector.eventNumber = self._particles[partPnt].event()
		self.numuDetector.pdgCode = self._particles[partPnt].pdgCode()
		self.numuDetector.x = self._particles[partPnt].x()
		self.numuDetector.y = self._particles[partPnt].y()
		self.numuDetector.z = self._particles[partPnt].z()
		self.numuDetector.s = self._particles[partPnt].s()
		self.numuDetector.px = self._particles[partPnt].p()[1][0]
		self.numuDetector.py = self._particles[partPnt].p()[1][1]
		self.numuDetector.pz = self._particles[partPnt].p()[1][2]
		self.numuDetector.t = self._particles[partPnt].t()
		self.numuDetector.eventWeight = self._particles[partPnt].weight()
		self.numuDetector.mass = self._particles[partPnt].mass()

# nue Detector
		partPnt = partPnt + 1
		self.nueDetector.runNumber = self._particles[partPnt].run()
		self.nueDetector.eventNumber = self._particles[partPnt].event()
		self.nueDetector.pdgCode = self._particles[partPnt].pdgCode()
		self.nueDetector.x = self._particles[partPnt].x()
		self.nueDetector.y = self._particles[partPnt].y()
		self.nueDetector.z = self._particles[partPnt].z()
		self.nueDetector.s = self._particles[partPnt].s()
		self.nueDetector.px = self._particles[partPnt].p()[1][0]
		self.nueDetector.py = self._particles[partPnt].p()[1][1]
		self.nueDetector.pz = self._particles[partPnt].p()[1][2]
		self.nueDetector.t = self._particles[partPnt].t()
		self.nueDetector.eventWeight = self._particles[partPnt].weight()
		self.nueDetector.mass = self._particles[partPnt].mass()

# write out
		self.evTree.Fill()

# empty the eventHistory once it has been written out
# don't zero, overwrite
#		for pnt in range(11):
#			self._particles[pnt] = None
		self.makeHistory()


#  Write the root structure out to the file
	def write(self):
		self.evTree.Write()

#  add particles to the history
	def addParticle(self, location, par):
		if location == "target":
			self._particles[0] = par
		elif location == "productionStraight":
			self._particles[1] = par
		elif location == "prodStraightEnd":
			self._particles[2] = par
		elif location == "pionDecay":
			self._particles[3] = par
		elif location == "muonProduction":
			self._particles[4] = par
		elif location == "piFlashNu":
			self._particles[5] = par
		elif location == "muonDecay":
			self._particles[6] = par
		elif location == "eProduction":
			self._particles[7] = par
		elif location == "numuProduction":
			self._particles[8] = par
		elif location == "nueProduction":
			self._particles[9] = par
		elif location == "numuDetector":
			self._particles[10] = par
		elif location == "nueDetector":
			self._particles[11] = par
		else:
			print ("unrecognised point on the history .. adding a particle", location)
			sys.exit(0)

#  find particles in the history
	def findParticle(self, location):
		if location == "target":
			return deepcopy(self._particles[0])
		elif location == "productionStraight":
			return deepcopy(self._particles[1])
		elif location == "prodStraightEnd":
			return deepcopy(self._particles[2])
		elif location == "pionDecay":
			return deepcopy(self._particles[3])
		elif location == "muonProduction":
			return deepcopy(self._particles[4])
		elif location == "piFlashNu":
			return deepcopy(self._particles[5])
		elif location == "muonDecay":
			return deepcopy(self._particles[6])
		elif location == "eProduction":
			return deepcopy(self._particles[7])
		elif location == "numuProduction":
			return deepcopy(self._particles[8])
		elif location == "nueProduction":
			return deepcopy(self._particles[9])
		elif location == "numuDetector":
			return deepcopy(self._particles[10])
		elif location == "nueDetector":
			return deepcopy(self._particles[11])

		else:
			print ("unrecognised point on the history .. finding a particle", location)
			sys.exit(0)

#  print out the eventHistory
	def display(self):

		print("Particle at target ", self._particles[0])
		print("Particle at start of production straight", self._particles[1])
		print("Particle at end of production straight", self._particles[2])
		print("Pion at decay point ", self._particles[3])
		print("muon at production point ", self._particles[4])
		print("neutrino from pion Flash ", self._particles[5])
		print("muon decay ", self._particles[6])
		print("electron Production ", self._particles[7])
		print("numu Production ", self._particles[8])
		print("nue Production ", self._particles[9])
		print("numu at Detector ", self._particles[10])
		print("nue at Detector ", self._particles[11])
