#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Class: eventHistory:
====================

    @file eventHistory.py

    @brief   creates a history of a nuStorm event

    @author  Paul Kyberd

    @version     1.0
    @date        16 September 2021

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
"struct atTarget {\
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
"struct prdStrght {\
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
		self._particles = [None]*11
		self.evTree = TTree('eventHistory', 'nuStorm Event')
		self._entryPnt = 0

	def __repr__(self):
		return "a history of a nuSTORM event"

	def __str__(self):
		return "eventHistory: run Number %i, outputFilename = %s" \
		         % \
               (self._runNum, self._outputFilename)


#	methods
#  Name and open a file for output
	def outFile(self, fileName ):
		self._outputFilename = fileName
		self._outTFile = TFile( self._outputFilename, 'RECREATE', 'event History' )

#  Name and open a file for input 
	def inFile(self, fileName ):
		self._inputFilename = fileName
		self._inTFile = TFile( self._inputFilename, 'READ', 'event History' )

		self._eHTree = self._inTFile.Get("eventHistory")
		self.atTarget = ROOT.atTarget()
		self._eHTree.SetBranchAddress('atTarget', self.atTarget)
		self.prodStrght = ROOT.prdStrght()
		self._eHTree.SetBranchAddress('prodStrght', self.prodStrght)
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
#			print ("at target x ", self.atTarget.x)
#			print ("production Straight x ", self.prodStrght.x)
#			print ("pion Decay x ", self.pionDecay.x)
#			print (self.pionDecay.runNumber)
#			print (self.pionDecay.eventNumber)

		self._eHTree.GetEntry(self._entryPnt)

		pTar = particle.particle(self.atTarget.runNumber, self.atTarget.eventNumber, self.atTarget.s, 
			self.atTarget.x, self.atTarget.y, self.atTarget.z, self.atTarget.px, self.atTarget.py, self.atTarget.pz, 
			self.atTarget.t, self.atTarget.eventWeight, self.atTarget.mass, self.atTarget.pdgCode)
		pPs = particle.particle(self.prodStrght.runNumber, self.prodStrght.eventNumber, self.prodStrght.s, 
			self.prodStrght.x, self.prodStrght.y, self.prodStrght.z, self.prodStrght.px, self.prodStrght.py, self.prodStrght.pz, 
			self.prodStrght.t, self.prodStrght.eventWeight, self.prodStrght.mass, self.prodStrght.pdgCode)
		pPiDcy = particle.particle(self.pionDecay.runNumber, self.pionDecay.eventNumber, self.pionDecay.s, 
			self.pionDecay.x, self.pionDecay.y, self.pionDecay.z, self.pionDecay.px, self.pionDecay.py, self.pionDecay.pz, 
			self.pionDecay.t, self.pionDecay.eventWeight, self.pionDecay.mass, self.pionDecay.pdgCode)
		pMuPrd = particle.particle(self.muonProduction.runNumber, self.muonProduction.eventNumber, self.muonProduction.s, 
			self.muonProduction.x, self.muonProduction.y, self.muonProduction.z, self.muonProduction.px, self.muonProduction.py, self.muonProduction.pz, 
			self.muonProduction.t, self.muonProduction.eventWeight, self.muonProduction.mass, self.muonProduction.pdgCode)
		pFlsNu = particle.particle(self.piFlashProd.runNumber, self.piFlashProd.eventNumber, self.piFlashProd.s, 
			self.piFlashProd.x, self.piFlashProd.y, self.piFlashProd.z, self.piFlashProd.px, self.piFlashProd.py, self.piFlashProd.pz, 
			self.piFlashProd.t, self.piFlashProd.eventWeight, self.piFlashProd.mass, self.piFlashProd.pdgCode)
		pMuDcy = particle.particle(self.muonDecay.runNumber, self.muonDecay.eventNumber, self.muonDecay.s, 
			self.muonDecay.x, self.muonDecay.y, self.muonDecay.z, self.muonDecay.px, self.muonDecay.py, self.muonDecay.pz, 
			self.muonDecay.t, self.muonDecay.eventWeight, self.muonDecay.mass, self.muonDecay.pdgCode)
		pEPrd = particle.particle(self.eProduction.runNumber, self.eProduction.eventNumber, self.eProduction.s, 
			self.eProduction.x, self.eProduction.y, self.eProduction.z, self.eProduction.px, self.eProduction.py, self.eProduction.pz, 
			self.eProduction.t, self.eProduction.eventWeight, self.eProduction.mass, self.eProduction.pdgCode)
		pNumuPrd = particle.particle(self.numuProduction.runNumber, self.numuProduction.eventNumber, self.numuProduction.s, 
			self.numuProduction.x, self.numuProduction.y, self.numuProduction.z, self.numuProduction.px, self.numuProduction.py, self.numuProduction.pz, 
			self.numuProduction.t, self.numuProduction.eventWeight, self.numuProduction.mass, self.numuProduction.pdgCode)
		pNuePrd = particle.particle(self.numuProduction.runNumber, self.numuProduction.eventNumber, self.nueProduction.s, 
			self.nueProduction.x, self.nueProduction.y, self.nueProduction.z, self.nueProduction.px, self.nueProduction.py, self.nueProduction.pz, 
			self.nueProduction.t, self.nueProduction.eventWeight, self.nueProduction.mass, self.nueProduction.pdgCode)
		pNumuDet = particle.particle(self.numuDetector.runNumber, self.numuDetector.eventNumber, self.numuDetector.s, 
			self.numuDetector.x, self.numuDetector.y, self.numuDetector.z, self.numuDetector.px, self.numuDetector.py, self.numuDetector.pz, 
			self.numuDetector.t, self.numuDetector.eventWeight, self.numuDetector.mass, self.numuDetector.pdgCode)
		pNueDet = particle.particle(self.nueDetector.runNumber, self.nueDetector.eventNumber, self.nueDetector.s, 
			self.nueDetector.x, self.nueDetector.y, self.nueDetector.z, self.nueDetector.px, self.nueDetector.py, self.nueDetector.pz, 
			self.nueDetector.t, self.nueDetector.eventWeight, self.nueDetector.mass, self.nueDetector.pdgCode)

		self.addParticle('target', pTar)
		self.addParticle('productionStraight', pPs)
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
		self.atTarget = ROOT.atTarget()
		self.evTree.Branch('atTarget', self.atTarget, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
		self.prodStrght = ROOT.prdStrght()
		self.evTree.Branch('prodStrght', self.prodStrght, 'run/I:event:pdgCode:x/F:y:z:s:px:py:pz:t:eventWeight:mass')
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

#  Close the output file
	def outFileClose(self ):
		self._outTFile.Close()

	def fill(self):

# target
		self.atTarget.runNumber = self._particles[0].run()
		self.atTarget.eventNumber = self._particles[0].event()
		self.atTarget.pdgCode = self._particles[0].pdgCode()
		self.atTarget.x = self._particles[0].x()
		self.atTarget.y = self._particles[0].y()
		self.atTarget.z = self._particles[0].z()
		self.atTarget.s = self._particles[0].s()
		self.atTarget.px = self._particles[0].p()[1][0]
		self.atTarget.py = self._particles[0].p()[1][1]
		self.atTarget.pz = self._particles[0].p()[1][2]
		self.atTarget.t = self._particles[0].t()
		self.atTarget.eventWeight = self._particles[0].weight()
		self.atTarget.mass = self._particles[0].mass()

# production straight
		self.prodStrght.runNumber = self._particles[1].run()
		self.prodStrght.eventNumber = self._particles[1].event()
		self.prodStrght.pdgCode = self._particles[1].pdgCode()
		self.prodStrght.x = self._particles[1].x()
		self.prodStrght.y = self._particles[1].y()
		self.prodStrght.z = self._particles[1].z()
		self.prodStrght.s = self._particles[1].s()
		self.prodStrght.px = self._particles[1].p()[1][0]
		self.prodStrght.py = self._particles[1].p()[1][1]
		self.prodStrght.pz = self._particles[1].p()[1][2]
		self.prodStrght.t = self._particles[1].t()
		self.prodStrght.eventWeight = self._particles[1].weight()
		self.prodStrght.mass = self._particles[1].mass()

# pion decay
		self.pionDecay.runNumber = self._particles[2].run()
		self.pionDecay.eventNumber = self._particles[2].event()
		self.pionDecay.pdgCode = self._particles[2].pdgCode()
		self.pionDecay.x = self._particles[2].x()
		self.pionDecay.y = self._particles[2].y()
		self.pionDecay.z = self._particles[2].z()
		self.pionDecay.s = self._particles[2].s()
		self.pionDecay.px = self._particles[2].p()[1][0]
		self.pionDecay.py = self._particles[2].p()[1][1]
		self.pionDecay.pz = self._particles[2].p()[1][2]
		self.pionDecay.t = self._particles[2].t()
		self.pionDecay.eventWeight = self._particles[2].weight()
		self.pionDecay.mass = self._particles[2].mass()

# muon production
		self.muonProduction.runNumber = self._particles[3].run()
		self.muonProduction.eventNumber = self._particles[3].event()
		self.muonProduction.pdgCode = self._particles[3].pdgCode()
		self.muonProduction.x = self._particles[3].x()
		self.muonProduction.y = self._particles[3].y()
		self.muonProduction.z = self._particles[3].z()
		self.muonProduction.s = self._particles[3].s()
		self.muonProduction.px = self._particles[3].p()[1][0]
		self.muonProduction.py = self._particles[3].p()[1][1]
		self.muonProduction.pz = self._particles[3].p()[1][2]
		self.muonProduction.t = self._particles[3].t()
		self.muonProduction.eventWeight = self._particles[3].weight()
		self.muonProduction.mass = self._particles[3].mass()

# neutrino flash production
		self.piFlashNu.runNumber = self._particles[4].run()
		self.piFlashNu.eventNumber = self._particles[4].event()
		self.piFlashNu.pdgCode = self._particles[4].pdgCode()
		self.piFlashNu.x = self._particles[4].x()
		self.piFlashNu.y = self._particles[4].y()
		self.piFlashNu.z = self._particles[4].z()
		self.piFlashNu.s = self._particles[4].s()
		self.piFlashNu.px = self._particles[4].p()[1][0]
		self.piFlashNu.py = self._particles[4].p()[1][1]
		self.piFlashNu.pz = self._particles[4].p()[1][2]
		self.piFlashNu.t = self._particles[4].t()
		self.piFlashNu.eventWeight = self._particles[4].weight()
		self.piFlashNu.mass = self._particles[4].mass()

# muon decay
		self.muonDecay.runNumber = self._particles[5].run()
		self.muonDecay.eventNumber = self._particles[5].event()
		self.muonDecay.pdgCode = self._particles[5].pdgCode()
		self.muonDecay.x = self._particles[5].x()
		self.muonDecay.y = self._particles[5].y()
		self.muonDecay.z = self._particles[5].z()
		self.muonDecay.s = self._particles[5].s()
		self.muonDecay.px = self._particles[5].p()[1][0]
		self.muonDecay.py = self._particles[5].p()[1][1]
		self.muonDecay.pz = self._particles[5].p()[1][2]
		self.muonDecay.t = self._particles[5].t()
		self.muonDecay.eventWeight = self._particles[5].weight()
		self.muonDecay.mass = self._particles[5].mass()

# e production
		self.eProduction.runNumber = self._particles[6].run()
		self.eProduction.eventNumber = self._particles[6].event()
		self.eProduction.pdgCode = self._particles[6].pdgCode()
		self.eProduction.x = self._particles[6].x()
		self.eProduction.y = self._particles[6].y()
		self.eProduction.z = self._particles[6].z()
		self.eProduction.s = self._particles[6].s()
		self.eProduction.px = self._particles[6].p()[1][0]
		self.eProduction.py = self._particles[6].p()[1][1]
		self.eProduction.pz = self._particles[6].p()[1][2]
		self.eProduction.t = self._particles[6].t()
		self.eProduction.eventWeight = self._particles[6].weight()
		self.eProduction.mass = self._particles[6].mass()

# numu production
		self.numuProduction.runNumber = self._particles[7].run()
		self.numuProduction.eventNumber = self._particles[7].event()
		self.numuProduction.pdgCode = self._particles[7].pdgCode()
		self.numuProduction.x = self._particles[7].x()
		self.numuProduction.y = self._particles[7].y()
		self.numuProduction.z = self._particles[7].z()
		self.numuProduction.s = self._particles[7].s()
		self.numuProduction.px = self._particles[7].p()[1][0]
		self.numuProduction.py = self._particles[7].p()[1][1]
		self.numuProduction.pz = self._particles[7].p()[1][2]
		self.numuProduction.t = self._particles[7].t()
		self.numuProduction.eventWeight = self._particles[7].weight()
		self.numuProduction.mass = self._particles[7].mass()

# nue production
		self.nueProduction.runNumber = self._particles[8].run()
		self.nueProduction.eventNumber = self._particles[8].event()
		self.nueProduction.pdgCode = self._particles[8].pdgCode()
		self.nueProduction.x = self._particles[8].x()
		self.nueProduction.y = self._particles[8].y()
		self.nueProduction.z = self._particles[8].z()
		self.nueProduction.s = self._particles[8].s()
		self.nueProduction.px = self._particles[8].p()[1][0]
		self.nueProduction.py = self._particles[8].p()[1][1]
		self.nueProduction.pz = self._particles[8].p()[1][2]
		self.nueProduction.t = self._particles[8].t()
		self.nueProduction.eventWeight = self._particles[8].weight()
		self.nueProduction.mass = self._particles[8].mass()

# numu Detector
		self.numuDetector.runNumber = self._particles[9].run()
		self.numuDetector.eventNumber = self._particles[9].event()
		self.numuDetector.pdgCode = self._particles[9].pdgCode()
		self.numuDetector.x = self._particles[9].x()
		self.numuDetector.y = self._particles[9].y()
		self.numuDetector.z = self._particles[9].z()
		self.numuDetector.s = self._particles[9].s()
		self.numuDetector.px = self._particles[9].p()[1][0]
		self.numuDetector.py = self._particles[9].p()[1][1]
		self.numuDetector.pz = self._particles[9].p()[1][2]
		self.numuDetector.t = self._particles[9].t()
		self.numuDetector.eventWeight = self._particles[9].weight()
		self.numuDetector.mass = self._particles[9].mass()

# nue Detector
		self.nueDetector.runNumber = self._particles[10].run()
		self.nueDetector.eventNumber = self._particles[10].event()
		self.nueDetector.pdgCode = self._particles[10].pdgCode()
		self.nueDetector.x = self._particles[10].x()
		self.nueDetector.y = self._particles[10].y()
		self.nueDetector.z = self._particles[10].z()
		self.nueDetector.s = self._particles[10].s()
		self.nueDetector.px = self._particles[10].p()[1][0]
		self.nueDetector.py = self._particles[10].p()[1][1]
		self.nueDetector.pz = self._particles[10].p()[1][2]
		self.nueDetector.t = self._particles[10].t()
		self.nueDetector.eventWeight = self._particles[10].weight()
		self.nueDetector.mass = self._particles[10].mass()

# write out
		self.evTree.Fill()

# empty the eventHistory once it has been written out
		for pnt in range(11):
			self._particles[pnt] = None

#  Write the root structure out to the file
	def write(self):
		self.evTree.Write()

#  add particles to the history
	def addParticle(self, location, par):
		if location == "target":
			self._particles[0] = par
		elif location == "productionStraight":
			self._particles[1] = par
		elif location == "pionDecay":
			self._particles[2] = par
		elif location == "muonProduction":
			self._particles[3] = par
		elif location == "piFlashNu":
			self._particles[4] = par
		elif location == "muonDecay":
			self._particles[5] = par
		elif location == "eProduction":
			self._particles[6] = par
		elif location == "numuProduction":
			self._particles[7] = par
		elif location == "nueProduction":
			self._particles[8] = par
		elif location == "numuDetector":
			self._particles[9] = par
		elif location == "nueDetector":
			self._particles[10] = par
		else:
			print ("unrecognised point on the history .. adding a particle", location)
			sys.exit(0)

#  find particles in the history
	def findParticle(self, location):
		if location == "target":
			return deepcopy(self._particles[0])
		elif location == "productionStraight":
			return deepcopy(self._particles[1])
		elif location == "pionDecay":
			return deepcopy(self._particles[2])
		elif location == "muonProduction":
			return deepcopy(self._particles[3])
		elif location == "piFlashNu":
			return deepcopy(self._particles[4])
		elif location == "muonDecay":
			return deepcopy(self._particles[5])
		elif location == "eProduction":
			return deepcopy(self._particles[6])
		elif location == "numuProduction":
			return deepcopy(self._particles[7])
		elif location == "nueProduction":
			return deepcopy(self._particles[8])
		elif location == "numuDetector":
			return deepcopy(self._particles[9])
		elif location == "nueDetector":
			return deepcopy(self._particles[10])

		else:
			print ("unrecognised point on the history .. finding a particle", location)
			sys.exit(0)

#  print out the eventHistory
	def display(self):

		print("Particle at target ", self._particles[0])
		print("Particle at start of production straight", self._particles[1])
		print("Pion at decay point ", self._particles[2])
		print("muon at production point ", self._particles[3])
		print("neutrino from pion Flash ", self._particles[4])
		print("muon decay ", self._particles[5])
		print("electron Production ", self._particles[6])
		print("numu Production ", self._particles[7])
		print("nue Production ", self._particles[8])
		print("numu at Detector ", self._particles[9])
		print("nue at Detector ", self._particles[10])

