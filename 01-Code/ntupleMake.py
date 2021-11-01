#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Class: ntupleMake:
==================

Outputs nuStorm monte-carlo production data as a root ntuple

Version 2.9                                    28/08/2021
Author: Ken Long
Update muon branch to include update to NeutrinoEventInstance class that
produces muon decays around the ring.

Version 2.8                                    26/05/2021
Author: Paul Kyberd
method to include information for pion flash into the flux section

Version 2.7                                    14/05/2021
Author: Paul Kyberd
Pion data goes in the electron location. Fixed

Version 2.6a                                    12/05/2021
Author: Paul Kyberd
Add a method to fill the ntuple from the pion flash data.
Don't change the version number; at present the pion data goes
in the electron location. Need to fix

Version 2.6									    26/02/2021
Author: Paul Kyberd
Add a branch for the downstream information


Version 2.5 									25/02/2021
Author: Paul Kyberd
Modify treeFill() method to output one event at a time

Version 2.4										24/02/2021
Author: Paul Kyberd
Include the time of the decay as part of the decay information.
This is the time in picoseconds since the start of the machine 
spill.


Version 2.3										17/02/2021
Author: Paul Kyberd
Make the run number a parameter passed in as part of __init__

Version 2.2										16/02/2021
Author: Paul Kyberd
Add run number to the data structure and include it in a new
branch of the tree, which is only written once at the start
of the run

Version 2.1										08/02/2021
Author: Paul Kyberd
Change the name of the root event object to reflect the version
of the format of the data output

Version 2.0										08/02/2021
Rewrite output as a TTree

Version 1.1										05/02/2021
Modify the variables output in the root file

Version 1.0										01/02/2021
Initial version

'''

# Generic Python imports
import sys, os

import array
# numpy
import numpy as np

# pyRoot imports
import ROOT
from ROOT import TFile, TNtuple, gROOT, TROOT, TTree, AddressOf, addressof

# nuStorm imports
import MuonConst as MuonConst
import PionConst as PionConst
import NeutrinoEventInstance as nuEvtInst

#  Need a structure so 
gROOT.ProcessLine(
"struct information {\
   Float_t          runNumber;\
   Float_t          Version;\
   Float_t			PWidth;\
   Float_t          PHeight;\
   Float_t          PZ;\
   Float_t			DBWidth;\
   Float_t          DBHeight;\
   Float_t          DBDepth;\
   Float_t          DBZ;\
   Char_t           emmitanceUnits[2];\
};" );
gROOT.ProcessLine(
"struct event_t {\
   Float_t         Epi;\
   Float_t         pPiX;\
   Float_t         pPiY;\
   Float_t         pPiZ;\
   Float_t         Emu;\
   Float_t         pMuX;\
   Float_t         pMuY;\
   Float_t         pMuZ;\
   Float_t         s;\
   Float_t         x;\
   Float_t         y;\
   Float_t         z;\
   Float_t         xp;\
   Float_t         yp;\
   Float_t         t;\
   Float_t         Ee;\
   Float_t         peX;\
   Float_t         peY;\
   Float_t         peZ;\
   Float_t         Enumu;\
   Float_t         pnumuX;\
   Float_t         pnumuY;\
   Float_t         pnumuZ;\
   Float_t         Enue;\
   Float_t         pnueX;\
   Float_t         pnueY;\
   Float_t         pnueZ;\
};" );

gROOT.ProcessLine(
"struct flux_t {\
   Float_t          nuEx;\
   Float_t          nuEy;\
   Float_t			nuEpx;\
   Float_t          nuEpy;\
   Float_t          nuEpz;\
   Float_t          nuEE;\
   Float_t			nuMux;\
   Float_t          nuMuy;\
   Float_t          nuMupx;\
   Float_t          nuMupy;\
   Float_t          nuMupz;\
   Float_t          nuMuE;\
};" );

muCnst = MuonConst.MuonConst()
piCnst = PionConst.PionConst()

class ntupleMake:
	Version = 2.6
	__Validated__ = False
	outputFilename="" 
	__muMass = muCnst.mass()/1000.
	__piMass = piCnst.mass()/1000.

# built-in methods
	def __init__(self, runNum, nuPrdStrt, outputFilename="nuStorm.root"):
# first a TTree with the run information - written once
		self.outputFilename = outputFilename
		self.outfile = TFile( self.outputFilename, 'RECREATE', 'ROOT file with an NTuple' )
# first a TTree with the run information - written once
		self.runInfo = ROOT.information()
		self.infoTree = TTree('runInfo', 'run information')
		self.infoTree.Branch( 'information', self.runInfo, 'runNumber/F:Version:FluxPlaneW:FluxPlaneH:PZ:DetectW:DetecH:DetecD:DetecZ:Emit')
#		self.infoTree.Branch( 'Flux Plane', AddressOf( self.runInfo, 'PWidth'), 'Width/F:Height:Z')
#		self.infoTree.Branch( 'Detector Box', AddressOf( self.runInfo, 'DBWidth'), 'Width/F:Height:Depth:Z')

# Parameters for the run
		self.runInfo.runNumber = runNum
		self.runInfo.Version = self.Version

# Parameters for the flux plane
		self.runInfo.PWidth = 10.0
		self.runInfo.PHeight = 10.0
		self.runInfo.PZ = nuPrdStrt.HallWallDist()
# Parameters for the detector box
		self.runInfo.DBWidth = nuPrdStrt.DetHlfWdth()*2.0
		self.runInfo.DBHeight = nuPrdStrt.DetHlfWdth()*2.0
		self.runInfo.DBDepth = nuPrdStrt.DetLngth()
		self.runInfo.DBZ = nuPrdStrt.Hall2Det()
		self.emittanceUnits = 'mm'

		self.infoTree.Fill()
		self.infoTree.Write()

# a TTree for the production data from every event
		self.event = ROOT.event_t()
		self.evTree = TTree('beam', 'nuStorm Event')
		self.evTree.Branch( 'pion', self.event, 'E/F:pX:pY:pZ' )
		self.evTree.Branch( 'muon', addressof( self.event, 'Emu'), 'E/F:pX:pY:pZ' )
		self.evTree.Branch( 'Decay', addressof( self.event, 's'), 's/F:x:y:z:xp:yp:t')
		self.evTree.Branch( 'Electron', addressof( self.event, 'Ee'), 'E/F:pX:pY:pZ')
		self.evTree.Branch( 'NuMu', addressof( self.event, 'Enumu'), 'E/F:pX:pY:pZ')
		self.evTree.Branch( 'Nue', addressof( self.event, 'Enue'), 'E/F:pX:pY:pZ')

# a TTree for the flux data
		self.flux = ROOT.flux_t()
		self.fluxTree = TTree('flux', 'nuStorm Event')
		self.fluxTree.Branch( 'NuE', addressof( self.flux, 'nuEx'), 'x/F:y:px:py:pz:E')
		self.fluxTree.Branch( 'NuMu', addressof( self.flux, 'nuMux'), 'x/F:y:px:py:pz:E')

	def __repr__(self):
		return "make an TTree"

	def __str__(self):
		return "ntuple Maker: outputFilename = %s" \
		         % \
               (self.outputFilename)

	def treeFill(self, nuEvt):

# tracespace coordinates of the decay point
		self.event.s = nuEvt.getTraceSpaceCoord()[0]
		self.event.x = nuEvt.getTraceSpaceCoord()[1]
		self.event.y = nuEvt.getTraceSpaceCoord()[2]
		self.event.z = nuEvt.getTraceSpaceCoord()[3]
		self.event.xp = nuEvt.getTraceSpaceCoord()[4]
		self.event.yp = nuEvt.getTraceSpaceCoord()[5]
		self.event.t = 0.0;
# get the 4 vector of the muon. In fact Emu ~ Pmu ~ PMuz - but calculate precisely
		pMu = nuEvt.getpmuGen()
		self.event.Emu  = np.sqrt(pMu*pMu + self.__muMass**2)
		pBeam = nuEvt.getPb()
		self.event.pMuX = pBeam[0]
		self.event.pMuY = pBeam[1]
		self.event.pMuZ = pBeam[2]
# get the 4 vector of the electron
		self.event.Ee  = nuEvt.gete4mmtm()[0]
		self.event.peX = nuEvt.gete4mmtm()[1][0]
		self.event.peY = nuEvt.gete4mmtm()[1][1]
		self.event.peZ = nuEvt.gete4mmtm()[1][2]
# get the 4 vector of the muon neutrino
		self.event.Enumu  = nuEvt.getnumu4mmtm()[0]
		self.event.pnumuX = nuEvt.getnumu4mmtm()[1][0]
		self.event.pnumuY = nuEvt.getnumu4mmtm()[1][1]
		self.event.pnumuZ = nuEvt.getnumu4mmtm()[1][2]
# get the 4 vector of the electron neutrino
		self.event.Enue  = nuEvt.getnue4mmtm()[0]
		self.event.pnueX = nuEvt.getnue4mmtm()[1][0]
		self.event.pnueY = nuEvt.getnue4mmtm()[1][1]
		self.event.pnueZ = nuEvt.getnue4mmtm()[1][2]

		self.evTree.Fill()

# fill the tree from a pion flash decau
	def pionTreeFill(self, piEvt):

# tracespace coordinates of the decay point
		self.event.s = piEvt.getTraceSpaceCoord()[0]
		self.event.x = piEvt.getTraceSpaceCoord()[1]
		self.event.y = piEvt.getTraceSpaceCoord()[2]
		self.event.z = piEvt.getTraceSpaceCoord()[3]
		self.event.xp = piEvt.getTraceSpaceCoord()[4]
		self.event.yp = piEvt.getTraceSpaceCoord()[5]
		self.event.t = 0.0;
# get the 4 vector of the decay muon.
		self.event.Emu  = piEvt.getmu4mmtm()[0]
		self.event.pMuX = piEvt.getmu4mmtm()[1][0]
		self.event.pMuY = piEvt.getmu4mmtm()[1][1]
		self.event.pMuZ = piEvt.getmu4mmtm()[1][2]
# get the 4 vector of the pion
		pPion = piEvt.getppiGen()
		self.event.Epi  = np.sqrt(pPion*pPion + self.__piMass**2)
		self.event.pPiX = pPion*self.event.xp
		self.event.pPiY = pPion*self.event.yp
		self.event.pPiZ = np.sqrt(pPion*pPion - self.event.peX**2 - self.event.peY**2)
# get the 4 vector of the muon neutrino
		self.event.Enumu  = piEvt.getnumu4mmtm()[0]
		self.event.pnumuX = piEvt.getnumu4mmtm()[1][0]
		self.event.pnumuY = piEvt.getnumu4mmtm()[1][1]
		self.event.pnumuZ = piEvt.getnumu4mmtm()[1][2]

		self.evTree.Fill()

	def fluxFill(self, hitE, hitMu):
		
		self.flux.nuEx = hitE[0]
		self.flux.nuEy = hitE[1]
		self.flux.nuEpx = hitE[5]
		self.flux.nuEpy = hitE[6]
		self.flux.nuEpz = hitE[7]
		self.flux.nuEE = hitE[8]
		self.flux.nuMux = hitMu[0]
		self.flux.nuMuy = hitMu[1]
		self.flux.nuMupx = hitMu[5]
		self.flux.nuMupy = hitMu[6]
		self.flux.nuMupz = hitMu[7]
		self.flux.nuMuE = hitMu[8]

		self.fluxTree.Fill()
#
#  Fill the flux branch from pion flash run - just remove the electron code
	def flashFluxFill(self, hitMu):

		self.flux.nuMux = hitMu[0]
		self.flux.nuMuy = hitMu[1]
		self.flux.nuMupx = hitMu[5]
		self.flux.nuMupy = hitMu[6]
		self.flux.nuMupz = hitMu[7]
		self.flux.nuMuE = hitMu[8]

		self.fluxTree.Fill()


	def output(self):
#	Write out the data
		self.evTree.Write()
		self.fluxTree.Write()
		return
		
		
	def closeFile(self):
		self.evTree.Print()
		self.evTree.Write()
		self.fluxTree.Print()
		self.fluxTree.Write()
		return

	def initNtuple(self, labels):
		self.ntuple  = TNtuple( 'ntuple', ' nuStorm source', ':'.join( labels ) )
		if self.__Validated__: print ("ntuple is ", self.ntuple)
		return


# Test to run if file is invoked as the main programme
if __name__ == "__main__" : 


	print(" running as main - run tests")

##! Generate events:
	nuEI = []
	for i in range(20):
		nuEI.append(nuEvtInst.NeutrinoEventInstance(5.))

	nt = ntupleMake("test.root")
	nt.treeFill(nuEI)	

	print(" tests finished - check file test.root")
