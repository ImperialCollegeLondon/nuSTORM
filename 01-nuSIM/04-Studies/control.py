#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class control:
==============

  Description of a pion

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
 1.0: 07Jan22: read the conditions for a particular run of the software
@author: PaulKyberd

 1.1: 06Jun22: include central stored muon energy and the ability to set a static runNumber
@author: MarvinPfaff

 1.2: 24May23: Include print flags to turn on and off printing of specific routines, but the code
               sets the flag to false if the corresponding flag is not defined in the dictionary file.
               This means that the code is backwards compatable and print flags do not have to be present
               in the dict file

 1.3: 03Sept23: Add a flag to set the x position of the detector centre

 1.4: 22Oct24:  Add a method to return the generation type

 1.5: 02Dec24:  Add a method for KDE generation to provide a directory to store the KDE pions



@author: PaulKyberd
"""

import math, sys
import os
import json
import numpy as np
from copy import deepcopy
import traceSpace
import MuonConst as mC
import PionConst as piC


class control:
    __Debug  = False

#--------  "Built-in methods":

    def __init__(self, controlFile):

        with open(controlFile) as controlFile:
            self._controlInfo = json.load(controlFile)

        self._var = 2.3
        self._runNumber = 0
# fill variables
        return

    def __repr__(self):
        return "control class"

    def __str__(self):
        return "study = %g,  mass = %g" % \
                  (self._var, self._var)

# Process decays in the transfer line
    def tlFlag(self):
        return (self._controlInfo["flags"]["tlFlag"] == "True")

# Process decays in the production straight
    def psFlag(self):
        return (self._controlInfo["flags"]["psFlag"] == "True")

# Process decays beyond the production straight
    def lstFlag(self):
        return (self._controlInfo["flags"]["lstFlag"] == "True")

# Process muon decays
    def muDcyFlag(self):
        return (self._controlInfo["flags"]["muDcyFlag"] == "True")

# Process muons outside the acceptance but in the production straight
    def PSMuons(self):
        return (self._controlInfo["flags"]["PSMuons"] == "True")
# Process muons that decay in the ring - after the overlap with the pions in the PS
    def ringMuons(self):
        return (self._controlInfo["flags"]["ringMuons"] == "True")
# Generate pencil beam at target
    def pencilBeam(self):
        return (self._controlInfo["flags"]["pencilBeam"] == "True")
# All pions start at t=0
    def tEqualsZero(self):
        return (self._controlInfo["flags"]["tEqualsZero"] == "True")
# Track the flash neutrinos to the detector
    def flashAtDetector(self):
        return (self._controlInfo["flags"]["flashAtDetector"] == "True")
# Get pion momentum distribution at target from histogram input
    def pDistInput(self):
        return (self._controlInfo["flags"]["pDistInput"] == "True")
# Get pion phase space distribution at target from histogram input
    def psDistInput(self):
        return (self._controlInfo["flags"]["psDistInput"] == "True")

# Get the genPion generation type
    def genType(self):
        return self._controlInfo["flags"]["genType"]

# Add possibility to set static runNumber
    def setRunNumber(self, runNum):
        self._runNumber = runNum
        print("========  Control Instance  ========")
        print("RunNumber set manually to ",self._runNumber)
        print("====================================")
        return True

# run number
    def runNumber(self, inc=False):

        if self._runNumber == 0:
    # run number is read from a file, it is in the studies directory a sub directory given by the study name and
    # a fileName given by the runNumber key word in the dictionary
    #    rNFile = "102-studies/" + self._controlInfo['study'] + "/" +self._controlInfo['runNumber']
            sDir =os.environ['StudyDir']
#        rNFile = "102-studies/" + self._controlInfo['study'] + "/" +self._controlInfo['runNumber']
            rNFile = sDir  + "/" +self._controlInfo['runNumber']
            rN = open(rNFile, "r")
            runNumber = int(rN.readline())
            rN.close()
            if (inc):
                runNumber = runNumber + 1
                rN = open(rNFile, "w")
                rN.write(str(runNumber))
                rN.close()
        else:
            runNumber = self._runNumber

        return runNumber

# number of events
    def nEvents(self):
        return self._controlInfo["nEvents"]

# x position of the centre of the detector in metres
    def detXPosition(self):
        try:
            xPos = self._controlInfo["detectorXPosition"]
            print(xPos)
        except KeyError as ke:
            print('Key Not Found in control dicitionary: providing zero as default', ke)
            xPos = 0.0
        return xPos

# print limit
    def printLimit(self):
        return self._controlInfo["printLimit"]

# Pion energy
    def PPi(self):
        return self._controlInfo["PPi"]

# Central muon energy stored in the ring
    def PMu(self):
        return self._controlInfo["PMu"]

#logFile name
    def logFile(self):
        sDir =os.environ['StudyDir'] + "/" + os.environ["StudyName"]
        return sDir + "/" + "/"+ self._controlInfo["files"]["logFile"] + str(self.runNumber()) + ".log"

#   Directory to put the fluka file in
    def flukaFileDir(self):
        return self._controlInfo["files"]["flukaFileDir"]


#plots dictionary file name
    def plotsDict(self):
        return self._controlInfo["files"]["plotsDict"]

#plots study name
    def studyName(self):
        return self._controlInfo["study"]

#run description study name
    def description(self):
        return self._controlInfo["description"]

#   Print out flags - main routine
    def mainPrnt(self):
        try:
            flag = self._controlInfo["print"]["main"]
        except KeyError as ke:
            flag = False
        print("in control mainPrnt flag set to ", flag)
        return flag

#   neutrino event instance
    def nuEvtInstPrnt(self):
        try:
            flag = self._controlInfo["print"]["nuEvtInst"]
        except KeyError as ke:
            flag = "False"
        return flag


