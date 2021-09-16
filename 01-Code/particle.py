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
  _PDG      : PDG code - so we track particle v anti-particle and any new stuff
    
  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates decay instance
      __repr__ : One liner with call.
      __str__  : Dump of values of decay

  Get/set methods:
  
  General methods:

Created on Tue 29Aug21;21:35: Version history:
----------------------------------------------
 1.0: 29Aug21: A straight copy of the pion class

@author: PaulKyberd
"""

import math
import numpy as np
from copy import deepcopy
import traceSpace
import abc
    
class particle:
    __Debug  = False
    
#--------  "Built-in methods":
    def __init__(self, runNum, eventNum, s, x, y, z, px, py, pz, t, eventWeight, mass, PDG):

        E = math.sqrt(px*px + py*py + pz*pz + mass*mass)
        p = np.array([px, py, pz])
        self._p = np.array([E, np.array([px, py, pz])],dtype=object)
        self._TrcSpc = traceSpace.traceSpace(s, x, y, z, px/pz, py/pz)
        self._t = t
        self._eventWeight = eventWeight
        self._mass = mass
        self._runNum = runNum
        self._eventNum = eventNum
        self._PDG = PDG
        return

    def __repr__(self):
        return "particle(x, y, z, s, px, py, pz, t, weight, mass)"

    def __str__(self):
        return "particle: E (GeV) = %g, p (GeV) = (%g  %g  %g),  mass = %g, t = %g, s = %g,\n x = %g, y = %g, z = %g, x' = %g, y' = %g, weight = %g " % \
                  (self._p[0], self._p[1][0],  self._p[1][1], self._p[1][2],self._mass, self._t, self._TrcSpc.s(), self._TrcSpc.x(), 
                    self._TrcSpc.y(), self._TrcSpc.z(), self._TrcSpc.xp() ,self._TrcSpc.yp(), self._eventWeight  )

# Override the __eq__ method
    def __eq__(self, comp):
        if isinstance(comp, self.__class__):
            if (self._p[0] == comp._p[0]) and (self._p[1][0] == comp._p[1][0]) and (self._p[1][1] == comp._p[1][1]) and \
               (self._p[1][2] == comp._p[1][2]) and (self._mass == comp._mass) and (self._t == comp._t)  and \
               (self._TrcSpc == comp._TrcSpc) and (self._eventWeight == comp._eventWeight) and (self._PDG == comp._PDG) and \
               (self._runNum == comp._runNum) and (self._eventNum == comp._eventNum):
                return True
            else:
                return False
        else:
            return False 

# Override the __ne__ method
    def __ne__(self, comp):
        return not self.__eq__(comp)

    
#--------  get/set methods:
    def run(self):
        return deepcopy(self._runNum)

    def event(self):
        return deepcopy(self._eventNum)

    def p(self):
        return deepcopy(self._p)

    def traceSpace(self):
        return deepcopy(self._TrcSpc)

    def t(self):
        return deepcopy(self._t)

    def weight(self):
        return deepcopy(self._eventWeight)

    def mass(self):
        return deepcopy(self._mass)

    def s(self):
        return deepcopy(self._TrcSpc.s())

    def x(self):
        return deepcopy(self._TrcSpc.x())

    def y(self):
        return deepcopy(self._TrcSpc.y())

    def z(self):
        return deepcopy(self._TrcSpc.z())

    def xp(self):
        return deepcopy(self._TrcSpc.xp())

    def yp(self):
        return deepcopy(self._TrcSpc.yp())
