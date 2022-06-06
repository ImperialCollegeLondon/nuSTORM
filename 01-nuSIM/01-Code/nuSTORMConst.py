#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class nuSTORMConst:
================

  Defines the constants that determine the properties of the nuSTORM design.
  The values represent current design values. Implemented as a singleton.

  Class attributes:
  -----------------
  __instance : Set on creation of first (and only) instance.

  Instance attributes:
  --------------------
  None; all constants returned by "get" methods.

  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __new__ : Creates singleton class and prints version,
                reference, and values of constants used.
      __repr__: One liner with call.
      __str__ : Dump of constants

  Get/set methods:
      CdVrsn()  : Returns code version number.
      PDGref()  : Returns reference to version of PDG used for
                  constants.
    pdgCode() : Returns the (absolute) value of the PDG code - 211
      mass()    : Mass (MeV)
      lifetime(): Lifetime (s)
      SoL()     : Speed of light (m/s)
      print()   : Dumps parameters

Created on Thu 26Jan2022; Version history:
----------------------------------------------
 1.0: 26Jan2022: First implementation - based on PionConst.py

@author: MarvinPfaff
"""

import os
import pandas as pnds
from copy import deepcopy

class nuSTORMConst(object):
    __instance = None

#--------  "Built-in methods":
    def __new__(cls):
        if cls.__instance is None:
            #print('nuSTRMCnst.__new__: creating the nuSTRMCnst object')
            cls.__instance = super(nuSTORMConst, cls).__new__(cls)
            cls._nuSIMPATH          = os.getenv('nuSIMPATH')
            cls._TrfLineFilename    = os.path.join(cls._nuSIMPATH, '11-Parameters/nuSTORM-TrfLineCmplx-Params-v1.0.csv')
            cls._PrdStrghtFilename  = os.path.join(cls._nuSIMPATH, '11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
            cls._TrfLineParams      = cls.GetParams(cls._TrfLineFilename)
            cls._TrfLineCmplxLen    = cls._TrfLineParams.iat[0,2]
            cls._TrfLineCmplxAng    = cls._TrfLineParams.iat[1,2]
            ##cls._piAcc              = cls._TrfLineParams.iat[2,2] / 100.
            ##cls._epsilon            = cls._TrfLineParams.iat[3,2]
            ##cls._beta               = cls._TrfLineParams.iat[4,2]
            cls._delT0              = cls._TrfLineParams.iat[5,2]
            cls._delT1              = cls._TrfLineParams.iat[6,2]
            cls._delT2              = cls._TrfLineParams.iat[7,2]
            cls._PrdStrghtParams    = cls.GetParams(cls._PrdStrghtFilename)
            cls._Circumference      = cls._PrdStrghtParams.iat[0,2]
            cls._ProdStrghtLen      = cls._PrdStrghtParams.iat[1,2]
            cls._piAcc              = cls._PrdStrghtParams.iat[2,2] / 100.
            cls._muAcc              = cls._PrdStrghtParams.iat[3,2] / 100.
            cls._epsilon            = cls._PrdStrghtParams.iat[4,2]
            cls._beta               = cls._PrdStrghtParams.iat[5,2]
            cls._HallWallDist       = cls._PrdStrghtParams.iat[6,2]
            cls._DetHlfWdth         = cls._PrdStrghtParams.iat[7,2]
            cls._DetLngth           = cls._PrdStrghtParams.iat[8,2]
            cls._Hall2Det           = cls._PrdStrghtParams.iat[9,2]
            cls._ArcLen             = cls._PrdStrghtParams.iat[10,2]

# Only constants; print values that will be used:
        #print("nuSTRMCnst: version:", cls.CdVrsn(cls))
        #print("nuSTRMCnst: PDG reference:", cls.PDGref(cls))
        #print("nuSTRMCnst: mass (MeV):", cls.mass(cls))
        #print("nuSTRMCnst: lifetime (s):", cls.lifetime(cls))
        #print("nuSTRMCnst: speed of light:", cls.SoL(cls))

        return cls.__instance

    def __repr__(self):
        return "nuSTORMConst()"

    def __str__(self):
        return "nuSTORMConst: version=%g, length of transfer line complex=%g m, \n " \
               "circumference=%g m, length of production straight=%g, Length of Arc=%g m, \n" \
               "Production straight to hall wall=%g m, detector halfwidth=%g, detector length=%g, hall wall to detector=%g m, \n" \
               "pi momentum acceptance=%g%%, mu momentum acceptance=%g%%, transverse acceptance=%g pi mm rad, beta=%g mm, \n" \
               "bunch length=%g ns, bunch spacing=%g ns, extraction length=%g us. \n" % \
               (self.CdVrsn(), self.TrfLineCmplxLen(), self.Circumference(), self.ProdStrghtLen(), self.ArcLen(), self.HallWallDist(), self.DetHlfWdth(), self.DetLngth(), self.Hall2Det(), \
                self.piAcc(), self.muAcc(), self.epsilon(), self.beta(), self.delT0(), self.delT1(), self.delT2() )

#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)

    def CdVrsn(self):
        return 1.0

#Transfer Line Parameters
    def TrfLineCmplxLen(self):
        return deepcopy(self._TrfLineCmplxLen)

    def TrfLineCmplxAng(self):
        return deepcopy(self._TrfLineCmplxAng)

    #def piAcc(self):
    #    return deepcopy(self._piAcc)

    #def epsilon(self):
    #    return deepcopy(self._epsilon)

    #def beta(self):
    #    return deepcopy(self._beta)

    def delT0(self):
        return deepcopy(self._delT0)

    def delT1(self):
        return deepcopy(self._delT1)

    def delT2(self):
        return deepcopy(self._delT2)

#Production Straight Parameters
    def Circumference(self):
        return deepcopy(self._Circumference)

    def ProdStrghtLen(self):
        return deepcopy(self._ProdStrghtLen)

    def ArcLen(self):
        return deepcopy(self._ArcLen)

    def piAcc(self):
        return deepcopy(self._piAcc)

    def muAcc(self):
        return deepcopy(self._muAcc)

    def epsilon(self):
        return deepcopy(self._epsilon)

    def beta(self):
        return deepcopy(self._beta)

    def HallWallDist(self):
        return deepcopy(self._HallWallDist)

    def DetHlfWdth(self):
        return deepcopy(self._DetHlfWdth)

    def DetLngth(self):
        return deepcopy(self._DetLngth)

    def Hall2Det(self):
        return deepcopy(self._Hall2Det)

#--------  I/o methods:
    def GetParams(_filename):
        Params = pnds.read_csv(_filename)
        return Params

#--------  Utilities:
    def printPrdStrghtParams(self):
        print(self._PrdStrghtParams)

    def printTrfLineParams(self):
        print(self._TrfLineParams)

    def print(self):
        print("nuSTRMCnst: Version:", self.CdVrsn())
        print("nuSTRMCnst: Length of transfer line complex (m):", self.TrfLineCmplxLen())
        print("nuSTRMCnst: Circumference of nuSTORM ring (m):", self.Circumference())
        print("nuSTRMCnst: Length of production straight (m):", self.ProdStrghtLen())
        print("nuSTRMCnst: Length of arc (m):", self.ArcLen())
        print("nuSTRMCnst: Distsance from production straight to hall wall (m):", self.HallWallDist())
        print("nuSTRMCnst: Detector halfwidth (m):", self.DetHlfWdth())
        print("nuSTRMCnst: Detector length (m):", self.DetLngth())
        print("nuSTRMCnst: Distance from hall wall to detector (m):", self.Hall2Det())
        print("nuSTRMCnst: Pion momentum acceptance:", self.piAcc())
        print("nuSTRMCnst: Muon momentum acceptance:", self.muAcc())
        print("nuSTRMCnst: Transverse acceptance (pi mm rad):", self.epsilon())
        print("nuSTRMCnst: Beta (mm):", self.beta())
        print("nuSTRMCnst: Bunch length (ns):", self.delT0())
        print("nuSTRMCnst: Bunch spacing (ns):", self.delT1())
        print("nuSTRMCnst: Length of extraction (us):", self.delT2())
