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
 2.0: 24Sep21: Unify all the particles and just distinguish with a pdg code
 1.1: 21Sep21: Add a constructor with the momentum vector
 1.0: 29Aug21: A straight copy of the pion class

@author: PaulKyberd
"""

import math, sys
import numpy as np
from copy import deepcopy
import traceSpace
import MuonConst as mC
import PionConst as piC

muCnst  = mC.MuonConst()
piCnst  = piC.PionConst()
    
class particle:
    __Debug  = False
    
#--------  "Built-in methods":
#  Constructor with elementary variables
#    def __init__(self, runNum, eventNum, s, x, y, z, px, py, pz, t, eventWeight, particleType):
    def __init__(self, *args):

        if (len(args) == 12):
# Uunpack args
            runNum = args[0]
            eventNum = args[1]
            s = args[2]
            x = args[3]
            y = args[4]
            z = args[5]
            px =args[6]
            py = args[7]
            pz = args[8]
            t = args[9]
            eventWeight = args[10]
            particleType = args[11]

# fill variables
            p = np.array([px, py, pz])
            self._TrcSpc = traceSpace.traceSpace(s, x, y, z, px/pz, py/pz)
            self._t = t
            self._eventWeight = eventWeight
            self._runNum = runNum
            self._eventNum = eventNum
        elif (len(args) == 10):
# Uunpack args
            runNum = args[0]
            eventNum = args[1]
            s = args[2]
            x = args[3]
            y = args[4]
            z = args[5]
            p =args[6]
            t = args[7]
            eventWeight = args[8]
            particleType = args[9]
# fill variables
            self._p = p
            self._TrcSpc = traceSpace.traceSpace(s, x, y, z, p[1][0]/p[1][2], p[1][1]/p[1][2])
            self._t = t
            self._eventWeight = eventWeight
            self._runNum = runNum
            self._eventNum = eventNum
        else:
            sys.exit("Number of arguments unrecognised " + str(len(args)))

# fill mass pdgCode and lifetime depending on the particle type
        if (particleType == "pi+"):
            mass = piCnst.mass()/1000.0
            lifetime = piCnst.lifetime()
            pdgCode = 211
        elif particleType == "pi-":
            mass = piCnst.mass()/1000.0
            lifetime = piCnst.lifetime()
            pdgCode = -211
        elif particleType == "mu-":
            mass = muCnst.mass()/1000.0
            lifetime = muCnst.lifetime()
            pdgCode = 13
        elif particleType == "mu+":
            mass = muCnst.mass()/1000.0
            lifetime = muCnst.lifetime()
            pdgCode = -13
        elif particleType == "mu-":
            mass = 0.105
            lifetime = 2.0*10E-06
            pdgCode = 13
        elif particleType == "e-":
            mass = 0.00511
            lifetime = math.inf
            pdgCode = 11
        elif particleType == "e+":
            mass = 0.00511
            lifetime = math.inf
            pdgCode = -11
        elif particleType == "nue":
            mass = 0.0
            lifetime = math.inf
            pdgCode = 12
        elif particleType == "nueBar":
            mass = 0.0
            lifetime = math.inf
            pdgCode = -12
        elif particleType == "numu":
            mass = 0.0
            lifetime = math.inf
            pdgCode = 14
        elif particleType == "numuBar":
            mass = 0.0
            lifetime = math.inf
            pdgCode = -14
        else:
            sys.exit("Unrecognised particle type")

        if (len(args) == 12):
            E = math.sqrt(px*px + py*py + pz*pz + mass*mass)
            self._p = np.array([E, p],dtype=object)

        self._mass = mass
        self._lifetime = lifetime
        self._PDG = pdgCode

        return

    def __repr__(self):
        return "particle(x, y, z, s, px, py, pz, t, weight, mass)"

    def __str__(self):
        return "particle: E (GeV) = %g, p (GeV) = (%g  %g  %g),  mass = %g, t = %g, s = %g,\n x = %g, y = %g, z = %g, x' = %g, y' = %g, \
eventweight = %g, run = %g, event = %g, PDG=%g" % \
                  (self._p[0], self._p[1][0],  self._p[1][1], self._p[1][2],self._mass, self._t, self._TrcSpc.s(), self._TrcSpc.x(), 
                    self._TrcSpc.y(), self._TrcSpc.z(), self._TrcSpc.xp() ,self._TrcSpc.yp(), self._eventWeight, self._runNum, self._eventNum, self._PDG )

# Override the __eq__ method
    def __eq__(self, comp):

        delta = 0.0001
        equals = True
        n = 0
        if isinstance(comp, self.__class__):
            if (abs(self._p[0] - comp._p[0]) > delta):
                equals = False
            if (abs(self._p[1][0] - comp._p[1][0]) > delta):
                equals = False
            if (abs(self._p[1][1] - comp._p[1][1]) > delta):
                equals = False
            if (abs(self._p[1][2] - comp._p[1][2]) > delta):
                equals = False
            if (abs(self._mass - comp._mass) > delta):
                equals = False
            if (abs(self._t - comp._t) > delta):
                equals = False
            if (abs(self._TrcSpc.s() - comp._TrcSpc.s()) > delta):
                equals = False
            if (abs(self._TrcSpc.x() - comp._TrcSpc.x()) > delta):
                equals = False
            if (abs(self._TrcSpc.y() - comp._TrcSpc.y()) > delta):
                equals = False
            if (abs(self._TrcSpc.z() - comp._TrcSpc.z()) > delta):
                equals = False
            if (abs(self._TrcSpc.xp() - comp._TrcSpc.xp()) > delta):
                equals = False
            if (abs(self._TrcSpc.yp() - comp._TrcSpc.yp()) > delta):
                equals = False
            if (abs(self._eventWeight - comp._eventWeight) >delta):
                equals = False
            if (abs(self._lifetime - comp._lifetime) >delta):
                equals = False
            if (self._runNum != comp._runNum):
                equals = False
            if (self._eventNum != comp._eventNum):
                equals = False
            if (self._PDG != comp._PDG):
                equals = False
            return equals
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

    def pdgCode(self):
        return deepcopy(self._PDG)

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
