
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class PionTimeDistribution:
================

  Generates time distributions for pion at target.

  Dependencies:
   - numpy, Simulation, math, deepcopy

  Instance attributes:
  --------------------
  _tb : bunch length
  _ts : bunch spacing, i.e. time between bunches
  _tx : duration of extraction

  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates time distribution instance:
                 Needed keyword arguments: t_bunch -- bunch length,
                 t_spacing -- bunch spacing, i.e. time between bunches,
                 t_extraction -- duration of one extraction
      __repr__ : One liner with call.
      __str__  : Dump of values of time distribution

  Get/set methods:
    getBunchLength      : Returns bunch length (s) of this time distribution instance
    getBunchSpacing     : Returns bunch spacing (s) of this time distribution instance
    getExtractionLength : Returns extraction duration (s) of this time distribution instance

  PionTimeDistribution methods:
    GenerateTime: Generates time of pion at target with respect to first pion.
                  Returns time (float). Units s

Created on Mon 22Nov21: Version history:
----------------------------------------------
 1.0: 22Nov21: First implementation

@author: MarvinPfaff
"""

import numpy as np
import math as math
from copy import deepcopy
import Simulation as Simu

class PionTimeDistribution:

#--------  "Built-in methods":
    def __init__(self,t_bunch,t_spacing,t_extraction):

        self._tb = t_bunch * 10**(-9)
        self._ts = t_spacing * 10**(-9)
        self._tx = t_extraction * 10**(-6)

        return

    def __repr__(self):
        return "PionTimeDistribution()"

    def __str__(self):
        return "PionTimeDistribution: bunch length: %g ns, \r\n     \
                        bunch spacing: %g ns, \r\n   \
                        extraction length: %g us \
               " % ( self._tb,              \
                     self._ts,              \
                     self._tx)

#--------  "Dynamic methods"; individual time
    def GenerateTime(self, **kwargs):
        check = False
        while not check:
            t = Simu.getRandom() * self._tx
            t_mod = math.fmod(t,self._tb+self._ts)
            if t_mod < self._tb:
                check = True
        return t

#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)
    def getBunchLength(self):
        return self._tb

    def getBunchSpacing(self):
        return self._ts

    def getExtractionLength(self):
        return self._tx
