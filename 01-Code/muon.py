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
  _p        : muon 4 momentum: (E, array(px, py, pz)), GeV
  _weight   : event Weight
  _mass     : pion mass
  _runNum   : run number
  _eventNum : event number
  
  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates decay instance
      __repr__ : One liner with call.
      __str__  : Dump of values of decay

  Get/set methods:
    run():          get the run number
    event():        get the event number
    p():            get the 4 vector
    traceSpace():   get the trace space vectr
    t():            get the time
    weight():       get the event weight
    s():            get s (trace space)
    x():            get x (trace space)
    y():            get y (trace space)
    z():            get z (trace space)
    xp():           get xp (trace space)
    yp():           get yp (trace space)
    
  General methods:

Created on Wed 01Sept 09:50: Version history:
----------------------------------------------
 1.0: 20Aug21: First implementation - attempting inheritance
 1.1: 14Sept21: Add run number 


@author: PaulKyberd
"""

from particle import particle
import MuonConst as muConst

muCnst = muConst.MuonConst()

class muon(particle):

# the generic particle has no mass defined so it must provided. It is done by
# getting the mass from the muonConstant class and using it to call the particle
# constructor.

    def __init__(self, runNum, eventNum, s, x, y, z, px, py, pz, t, weight):
        mass = muCnst.mass()/1000.
        muonPDG = 13
        particle.__init__(self, runNum, eventNum, s, x, y, z, px, py, pz, t, weight, mass,muonPDG)

    def __repr__(self):
        return "particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight, mass)"


