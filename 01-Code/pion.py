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
    
  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates decay instance
      __repr__ : One liner with call.
      __str__  : Dump of values of decay

  Get/set methods:
  
  General methods:

Created on Tue 20Aug12;35: Version history:
----------------------------------------------
 1.0: 20Aug21: First implementation
 1.1: 22Aug21: Use the traceSpace class

@author: PaulKyberd
"""

import math
import numpy as np
from copy import deepcopy
import PionConst as PionConst
import traceSpace

piCnst = PionConst.PionConst()
    
class pion:

    __pimass = piCnst.mass()/1000.
    __Debug  = False
    
#--------  "Built-in methods":
    def __init__(self, s, x, y, z, px, py, pz, t, weight):

        E = math.sqrt(px*px + py*py + pz*pz + self.__pimass*self.__pimass)
        p = np.array([px, py, pz])
        self._p = np.array([E, np.array([px, py, pz])],dtype=object)
        self._TrcSpc = traceSpace.traceSpace(s, x, y, z, px/pz, py/pz)
        self._t = t
        self._weight = weight
        return

    def __repr__(self):
        return "pion(x, y, z, s, px, py, pz, t, weight)"

    def __str__(self):
        return "pion: E (GeV) = %g, p (GeV) = (%g  %g  %g),  t = %g, s = %g, x = %g, y = %g, z = %g, x' = %g, y' = %g, weight = %g " % \
                  (self._p[0], self._p[1][0],  self._p[1][1], self._p[1][2],self._t, self._TrcSpc.s(), self._TrcSpc.x(), 
                    self._TrcSpc.y(), self._TrcSpc.z(), self._TrcSpc.xp() ,self._TrcSpc.yp(), self._weight  )
    
#--------  get/set methods:
    def getP(self):
        return deepcopy(self._p)

    def traceSpace(self):
        return deepcopy(self._TrcSpc)

    def t(self):
        return deepcopy(self._t)

    def weight(self):
        return deepcopy(self._weight)

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
