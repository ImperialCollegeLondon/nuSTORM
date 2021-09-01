#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class muon:
===========

  Description of a muon

  Class attributes:

  --np       : numpy class
      
  Instance attributes:
  --------------------
  _TrcSpcCrd: Trace space (s, x, y, z, x', y') in numpy array - what are s and z?
  _t        : time in nanoseconds
  _p        : Pion 4 momentum: (E, array(px, py, pz)), GeV
  _weight   : not sure how to define
    
  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates decay instance
      __repr__ : One liner with call.
      __str__  : Dump of values of decay

  Get/set methods:
  
  General methods:

Created on Wed 01Sept 09:50: Version history:
----------------------------------------------
 1.0: 20Aug21: First implementation - attempting inheritance


@author: PaulKyberd
"""

from particle import particle
import MuonConst as muConst

muCnst = muConst.MuonConst()

class muon(particle):

# the generic particle has no mass defined so it must provided. It is done by
# getting the mass from the muonConstant class and using it to call the particle
# constructor.

    def __init__(self, s, x, y, z, px, py, pz, t, weight):
        mass = muCnst.mass()/1000.
        particle.__init__(self, s, x, y, z, px, py, pz, t, weight, mass)

    def __repr__(self):
        return "particle(x, y, z, s, px, py, pz, t, weight)"


