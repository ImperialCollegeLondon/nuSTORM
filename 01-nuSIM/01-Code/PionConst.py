#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class PionConst:
================

  Defines the constants that determine the properties of the muon.
  The values are taken from the PDG. Implemented as a singleton

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
    pdgCode() : Returns the (absolute) value of the PDG code - 211
      mass()    : Mass (MeV)
      lifetime(): Lifetime (s)
      SoL()     : Speed of light (m/s)
      print()   : Dumps parameters

Created on Thu 25Mar2021; Version history:
----------------------------------------------
 1.0: 25Mar2021: First implementation - based on MuonConst.py

@author: kennethlong and Paul Kyberd
"""

class PionConst(object):
    __instance = None

#--------  "Built-in methods":
    def __new__(cls):
        if cls.__instance is None:
            #print('PionConst.__new__: creating the PionConst object')
            cls.__instance = super(PionConst, cls).__new__(cls)

# Only constants; print values that will be used:
        #print("PionConst: version:", cls.CdVrsn(cls))
        #print("PionConst: PDG reference:", cls.PDGref(cls))
        #print("PionConst: mass (MeV):", cls.mass(cls))
        #print("PionConst: lifetime (s):", cls.lifetime(cls))
        #print("PionConst: speed of light:", cls.SoL(cls))

        return cls.__instance

    def __repr__(self):
        return "PionConst()"

    def __str__(self):
        return "PionConst: version=%g, mass=%g MeV, lifetime=%g s, c=%g m/s, PDG=%s" % (self.CdVrsn(), self.mass(), self.lifetime(), self.SoL(), self.PDGref())

#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)

    def CdVrsn(self):
        return 1.0

    def PDGref(self):
        return "M. Tanabashi et al. (Particle Data Group), Phys. Rev. D 98, 030001 (2018)"

    def pdgCode(self):
        return 211

    def mass(self):
        return 139.57061

    def lifetime(self):
        return 2.6033E-8

    def SoL(self):
        return 299792458.

#--------  Utilities:
    def print(self):
        print("PionConst: version:", self.CdVrsn())
        print("PionConst: PDG reference:", self.PDGref())
        print("PionConst: PDG code", self.pdgCode())
        print("PionConst: mass (MeV):", self.mass())
        print("PionConst: lifetime (s):", self.lifetime())
        print("PionConst: speed of light (m/s):", self.SoL())
