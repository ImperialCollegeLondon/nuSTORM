#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class pion:
===========

  Description of a pion

  Class attributes:

  --np       : numpy class
      
  Instance attributes:
  --------------------
  _TrcSpcCrd: Trace space (s, x, y, z, x', y') in numpy array - what are s and z?
  _t        : time in nanoseconds
  _p        : Pion 4 momentum: (E, array(px, py, pz)), GeV
  _weight   : not sure how to define
  _mass     : pion mass
    
  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates decay instance
      __repr__ : One liner with call.
      __str__  : Dump of values of decay

  Get/set methods:
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

Created on Tue 20Aug12;35: Version history:
----------------------------------------------
 1.0: 20Aug21: First implementation
 1.1: 22Aug21: Use the traceSpace class
 1.2: 01Sep21: replace getP() by p() to return vector, for consistency
 2.0: 01Sep21: make the class a derived class of the particle class

@author: PaulKyberd
"""

from particle import particle
import PionConst as PionConst

piCnst = PionConst.PionConst()

class pion(particle):

# the generic particle has no mass defined so it must provided. It is done by
# getting the mass from the piConstant class and using it to call the particle
# constructor.

    def __init__(self, run, event, s, x, y, z, px, py, pz, t, weight):
        mass = piCnst.mass()/1000.
        particle.__init__(self, run, event, s, x, y, z, px, py, pz, t, weight, mass)

    def __repr__(self):
        return "particle(run, event, x, y, z, s, px, py, pz, t, weight)"
