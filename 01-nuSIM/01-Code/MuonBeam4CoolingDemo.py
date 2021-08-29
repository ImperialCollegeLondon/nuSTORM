#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class MuonBeam4CoolingDemo
==========================

  Defines the muon beam for cooling demo parameters and allows generation
  of muons in yhe beam.

  Reads .csv file containing necessary parameters.  At present the beam
  is sampled as independent Gaussians in x, y, p_x, p_y, t, and E.

  Class attributes:
  -----------------
  __instance : Set on creation of first (and only) instance.
      
  Instance attributes:
  --------------------
  _filename        = Filename, includin path, to cls file containing
                     muon beam for cooling demo specification
  _PrdStrghtParams - Pandas dataframe instance read from _filename
  _TransEmit       = Transverse emittance (mm)
  _TransBeta       = Transverse beta (mm)
  _TransSigma      = Transverse sigma, position (mm)
  _TransMmtm       = Transverse momentum sigma (MeV/c)
  _MmtmBar         = Mean momentum (MeV/c)
  _LongiEmit       = Longitudinal emittance (mm)
  _TimeSigma       = Time sigma (ns)
  _EnergySigma     = Energy sigma (ns)

  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __new__ : Creates singleton class and prints version and values
                of decay-straight parameters used.
      __repr__: One liner with call.
      __str__ : Dump of constants

  Get/set methods:
      CdVrsn()     : Returns code version number.

      FileName()   : Returns csv paramter file name with path
    
      TransEmit()  : Returns transverse emittance
    
      TransBeta()  : Returns transverse beta
    
      TransSigma() : Returns transverse position sigma
    
      TransMmtm()  : Returns transverse mmtm sigma
    
      MmtmBar()    : Returns mean momentum
    
      LongiEmit()  : Reeturns longitudinal emittance
    
      TimeSigma()  : Returns time sigma
    
      EnergySigma(): Returns energy sigma

  Generation, calcution methods and utilities:
      None so far

  
Created on Mon 17May21. Version history:
----------------------------------------
 1.0: 17May21: First implementation

@author: kennethlong
"""

import pandas as pnds
import numpy as np

class MuonBeam4CoolingDemo(object):
    __instance = None

#--------  "Built-in methods":
    def __new__(cls, filename="Dummy"):
        if cls.__instance is None:
            cls.__instance = super(MuonBeam4CoolingDemo, cls).__new__(cls)

            cls._filename        = filename
            cls._PrdStrghtParams = cls.GetClngDmoBmPrm(filename)

            cls._TransEmit   = cls._PrdStrghtParams.iat[0,2]
            cls._TransBeta   = cls._PrdStrghtParams.iat[1,2]
            cls._TransSigma  = cls._PrdStrghtParams.iat[2,2]
            cls._TransMmtm   = cls._PrdStrghtParams.iat[3,2]
            cls._MmtmBar     = cls._PrdStrghtParams.iat[4,2]
            cls._LongiEmit   = cls._PrdStrghtParams.iat[5,2]
            cls._TimeSigma   = cls._PrdStrghtParams.iat[6,2]
            cls._EnergySigma = cls._PrdStrghtParams.iat[7,2]

        return cls.__instance

    def __repr__(self):
        return "MuonBeam4CoolingDemo()"

    def __str__(self):
        return "MuonBeam4CoolingDemo: version=%g \n" \
               "    filename=%s \n" \
               "    TransEmit=%g mm," \
               " TransBeta=%g mm," \
               " TransSigma=%g mm," \
               " TransMmtm=%g MeV/c" \
               " MmtmBar=%g MeV/c," \
               " LongiEmit=%g mm, " \
               " TimeSigma=%g ns," \
               " EnergySigma=%g MeV" \
               % \
               (self.CdVrsn(), \
                self.FileName(), \
                self.TransEmit(), \
                self.TransBeta(), \
                self.TransSigma(), \
                self.TransMmtm(), \
                self.MmtmBar(), \
                self.LongiEmit(), \
                self.TimeSigma(), \
                self.EnergySigma() \
                )

#--------  I/o methods:
    def GetClngDmoBmPrm(_filename):
        Muons4CoolingDemo = pnds.read_csv(_filename)
        return Muons4CoolingDemo
    
#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)

    def CdVrsn(self):
        return 1.0

    def FileName(self):
        return self._filename
    
    def TransEmit(self):
        return self._TransEmit
    
    def TransBeta(self):
        return self._TransBeta
    
    def TransSigma(self):
        return self._TransSigma
    
    def TransMmtm(self):
        return self._TransMmtm
    
    def MmtmBar(self):
        return self._MmtmBar
    
    def LongiEmit(self):
        return self._LongiEmit
    
    def TimeSigma(self):
        return self._TimeSigma
    
    def EnergySigma(self):
        return self._EnergySigma
    
