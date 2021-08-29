#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class MuonConst:
================

  Defines the constants that determine the properties of the muon.
  The values are taken from the PDG.  Goal is to implement this class
  as a singleton class as an excersise and so that there is no
  ambiguity about which values are in use.

  Class attributes:
  -----------------
  __instance : Set on creation of first (and only) instance.
      
  Instance attributes:
  --------------------
  None; all constants returned by "get" methods.
    
  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __new__ : Creates singleton class and prints version, PDG
                reference, and values of constants used.
      __repr__: One liner with call.
      __str__ : Dump of constants

  Get/set methods:
      CdVrsn()  : Returns code version number.
      PDGref()  : Returns reference to version of PDG used for
                  constants.
      pdgCode() : Returns the (absolute) value of the PDG code - 13
      mass()    : Mass (MeV)
      lifetime(): Lifetime (s)
      Michel()  : Michel parameters -- list, 3
      SoL()     : Speed of light (m/s)
      print()   : Dumps parameters
  
Created on Thu 31Dec20;16:42: Version history:
----------------------------------------------
 1.1: 12Oct21: Add pdgcode (absolute value)
 1.0: 31Dec20: First implementation

@author: kennethlong
"""

class MuonConst(object):
    __instance = None

#--------  "Built-in methods":
    def __new__(cls):
        if cls.__instance is None:
            #print('MuonConst.__new__: creating the MuonConst object')
            cls.__instance = super(MuonConst, cls).__new__(cls)
            
# Only constants; print values that will be used:
        #print("MuonConst: version:", cls.CdVrsn(cls))
        #print("MuonConst: PDG reference:", cls.PDGref(cls))
        #print("MuonConst: mass (MeV):", cls.mass(cls))
        #print("MuonConst: lifetime (s):", cls.lifetime(cls))
        #print("MuonConst: SM Michel parameters [rho, eta, delta]:", cls.Michel(cls))
        #print("MuonConst: speed of light:", cls.SoL(cls))
            
        return cls.__instance

    def __repr__(self):
        return "MuonConst()"

    def __str__(self):
        return "MuonConst: version=%g, mass=%g MeV, lifetime=%g s, Michel parameters=[%g, %g, %g], c=%g m/s, PDG=%s" % (self.CdVrsn(), self.mass(), self.lifetime(), self.Michel()[0], self.Michel()[1], self.Michel()[2], self.SoL(), self.PDGref())

#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)

    def CdVrsn(self):
        return 1.0

    def PDGref(self):
        return "P.A. Zyla et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2020, 083C01 (2020)."

    def pdgCode(self):
        return 13

    def mass(self):
        return 105.6583745

    def lifetime(self):
        return 2.1969811E-6

    def Michel(self):
        return [0.75, 0.0, 0.75]

    def SoL(self):
        return 299792458.

#--------  Utilities:
    def print(self):
        print("MuonConst: version:", self.CdVrsn())
        print("MuonConst: PDG reference:", self.PDGref())
        print("MuonConst: PDG code", self.pdgCode())
        print("MuonConst: mass (MeV):", self.mass())
        print("MuonConst: lifetime (s):", self.lifetime())
        print("MuonConst: SM Michel parameters [rho, eta, delta]:", self.Michel())
        print("MuonConst: speed of light (m/s):", self.SoL())
        
