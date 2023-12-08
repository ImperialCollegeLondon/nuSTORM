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
#   Plot energy v mu position
        hTitle = "munu Energy v Radius"
        hBins1 = 100
        hLower1 = 0.0
        hUpper1 = 5.0
        hBins2 = 100
        hLower2 = 0.0
        hUpper2 = 20.0
        self.hEnumuPosition = self.anHm.book2(hTitle, hBins1, hLower1, hUpper1, hBins2, hLower2, hUpper2)

#   Plot energy v e position
        hTitle = "nue Energy v Radius"
        hBins1 = 100
        hLower1 = 0.0
        hUpper1 = 5.0
        hBins2 = 100
        hLower2 = 0.0
        hUpper2 = 20.0
        self.hEnuePosition = self.anHm.book2(hTitle, hBins1, hLower1, hUpper1, hBins2, hLower2, hUpper2)

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
        if ((abs(hitPosX) < 2.5) and (abs(hitPosY) < 2.5)):
            numuP = numuDet.p()
            ENumuSpectra = numuP[0]
            self.hNumuDet.Fill(ENumuSpectra)

#   Fill the scatter plot for all entries
        radiusMu = math.sqrt(hitPosX*hitPosX + hitPosY*hitPosY)
        self.hEnumuPosition.Fill(numuDet.p()[0],radiusMu)

        return

    def processNueDet(self, nueDet):

        hitPosX = nueDet.x() - self.detPos[0]
        hitPosY = nueDet.y() - self.detPos[1]
        if ((abs(hitPosX) < 2.5) and (abs(hitPosY) < 2.5)):
            nueP = nueDet.p()
            ENueSpectra = nueP[0]
            self.hNueDet.Fill(ENueSpectra)

#   Fill the scatter plot for all entries
        radiusE = math.sqrt(hitPosX*hitPosX + hitPosY*hitPosY)
        self.hEnuePosition.Fill(nueDet.p()[0],radiusE)

        return



    def conclude(self, fileName):

      self.anHm.histOutRoot(fileName)

      return

