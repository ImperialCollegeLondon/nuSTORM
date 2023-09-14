#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class analysis:
==============

  Analyse an event during the generation process

  Class attributes:



  Instance attributes:
  --------------------


  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates a instance of the control class
      __repr__ : One liner with call.
      __str__  : Dump the conditions

  Get/set methods:

  General methods:

Version history:
----------------
 1.0: 11Sept23: Get a simple plot working

@author: PaulKyberd


"""

import math
import histoManager
import particle

class analyse:

    Debug  = True

#--------  "Built-in methods":

    def __init__(self, detectorPosition):

      #  need to store the detector position
        self.detPos = detectorPosition
        if (analyse.Debug) : print (f"detPos is {self.detPos}")

      # book a histogram
        self.anHm = histoManager.histoManager()
        hBins = 360
        hLower = 0.0
        hUpper = 7.2
        hTitle = "Energy of muon neutrinos at the detector plane"
        self.hNumuDet = self.anHm.book(hTitle, hBins, hLower, hUpper)
        hTitle = "Energy of electron neutrinos at the detector plane"
        self.hNueDet = self.anHm.book(hTitle, hBins, hLower, hUpper)

        self._var = 1.0

        return

    def __repr__(self):
        return "analyse class"

    def __str__(self):
        return "version is = %g" % \
                  (self._var)

    def processNumuDet(self, numuDet):

# Calculate the position of the hit on the front face
#          numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, pxnu, pynu, pznu, tNumu, eW, "numu")

        hitPosX = numuDet.x() - self.detPos[0]
        hitPosY = numuDet.y() - self.detPos[1]
        print(" ++ analyse: testing numu - x,y is ", hitPosX, ",  ", hitPosY)
        if ((abs(hitPosX) < 2.5) and (abs(hitPosY) < 2.5)):
            numuP = numuDet.p()
            ENumuSpectra = numuP[0]
            print("       ++ analyse: filling numu spectra - energyis ", ENumuSpectra)
            self.hNumuDet.Fill(ENumuSpectra)
        else:
            print("       ++ analyse: Not filling nuMu")


        return

    def processNueDet(self, nueDet):

        hitPosX = nueDet.x() - self.detPos[0]
        hitPosY = nueDet.y() - self.detPos[1]
        print(" ++ analyse: testing nue - x,y is ", hitPosX, ",  ", hitPosY)
        if ((abs(hitPosX) < 2.5) and (abs(hitPosY) < 2.5)):
            nueP = nueDet.p()
            ENueSpectra = nueP[0]
            print("       ++ analyse: filling nue spectra - energyis ", ENueSpectra)
            self.hNueDet.Fill(ENueSpectra)
        else:
            print("       ++ analyse: Not filling nuE")


        return



    def conclude(self, fileName):

      self.anHm.histOutRoot(fileName)

      return

