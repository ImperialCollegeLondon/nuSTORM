#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Class: Plots:
=============

Plot beam parameters for numus from the pion flash

Version 1.0										07/07/2021
PK - Somewhere to put root histograms


'''

# Generic Python imports
import os, time
import math

# pyRoot imports
from ROOT import gRandom
from ROOT import TFile
import ROOT

class Plots:
	__version__ = 1.0
	__Validated__ = False
	
	# Initialise the plots
	def __init__(self):
		print ("Initialise plots ")

		# Define the histograms ... here we have an example which is the energy of the muon neutrino
		self.title = "E numu"
		name = self.title								# name must be unique .. put it equal to title
		bins = 100
		lowerX =  0.0 			# range plus and minus 100m
		upperX =  10
		self.eNumuPlt = ROOT.TH1D(name, self.title, bins, lowerX, upperX)
		
		# add plots here
		return
		
	def fill(self, hitMu):
		
		nuMuE = hitMu[8]
		
#
#  Pass an object which can supply values via a get method        
# Fill the histogram with the values
		self.eNumuPlt.Fill(nuMuE)
		
		return
		
	def histdo(self):
		print("Histdo")
		# output the histograms to a root file
		self.outfile = TFile( 'plots.root', 'RECREATE', 'ROOT file with Histograms' )
		self.eNumuPlt.Write()
		
		return

if __name__ == "__main__" : 

	print (" running as main - run tests on plots. Version ", Plots.__version__)

