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
"""

import math, sys
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

    def ringMuons(self):
        return (self._controlInfo["flags"]["ringMuons"] == "True")


# Track the flash neutrinos to the detector
    def flashAtDetector(self):
        return (self._controlInfo["flags"]["flashAtDetector"] == "True")

# run number 
    def runNumber(self, inc=False):

    # run number is read from a file, it is in the studies directory a sub directory given by the study name and
    # a fileName given by the runNumber key word in the dictionary
        rNFile = "101-studies/" + self._controlInfo['study'] + "/" +self._controlInfo['runNumber']
        rN = open(rNFile, "r")
        runNumber = int(rN.readline())
        rN.close()
        if (inc):
            runNumber = runNumber + 1
            rN = open(rNFile, "w")
            rN.write(str(runNumber))
            rN.close()
        return runNumber

# number of events
    def nEvents(self):
        return self._controlInfo["nEvents"]

# Pion energy
    def EPi(self):
        return self._controlInfo["EPi"]

#logFile name
    def logFile(self):
        return "101-studies/" + self._controlInfo['study'] + "/"+ self._controlInfo["files"]["logFile"] + str(self.runNumber()) + ".log"

#plots dictionary file name
    def plotsDict(self):
        return self._controlInfo["files"]["plotsDict"]

#plots study name
    def studyName(self):
        return self._controlInfo["study"]

#run description study name
    def description(self):
        return self._controlInfo["description"]
