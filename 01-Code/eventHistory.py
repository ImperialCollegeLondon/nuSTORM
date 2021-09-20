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

class eventHistory:
	Version = 1.0
	__Validated__ = False

# built in methods
	def __init__(self ):
		self._outputFilename = "null"
		self._runNum = -1
		self._particles = [None]*11

	def __repr__(self):
		return "a history of a nuSTORM event"

	def __str__(self):
		return "eventHistory: run Number %i, outputFilename = %s" \
		         % \
               (self._runNum, self._outputFilename)


#	methods
#  Name and open the output file
	def outFile(self, fileName ):
		self._outputFilename = fileName
		self._outTFile = TFile( self._outputFilename, 'RECREATE', 'event History' )
#  Close the output file
	def outFileClose(self ):
		self._outTFile.Close()

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

