#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class mu4ClDmo
==============

  Generates an instance of a muon in the beam sent to the cooling demo.

  Assumes beam parameters already available through the MuonBeam4ClDmo 
  class

  Class attributes:
  -----------------
  __Debug : Controls print out
      
  Instance attributes:
  --------------------
  - _Rmu: 4-vector [t, np.array(x,y,z)]; ns, mm
  - _Pmu: 4-vector [E, np.array(px, py, pz)], MeV, MeV/c

  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __new__ : Creates instance
      __repr__: One liner with call.
      __str__ : Dump of constants

  Get/set methods:
      CdVrsn()     : Returns code version number.
"""

import MuonBeam4CoolingDemo
import MuonConst as muConst
import numpy as np

mc         = muConst.MuonConst()
muBm4ClDmo = MuonBeam4CoolingDemo.MuonBeam4CoolingDemo()

#--------  Class to get beam:
class mu4ClDmo:
    __Debug = False

    def __init__(self):
        if mu4ClDmo.__Debug:
            print(" mu4ClDmo.__init__:")

        Rmu = [0., np.array([0., 0., 0.])]
        Pmu = [0., np.array([0., 0., 0.])]

        Rmu[0]    = np.random.normal(0., muBm4ClDmo.TimeSigma(), None)
        Rmu[1][0] = np.random.normal(0., muBm4ClDmo.TransSigma(), None)
        Rmu[1][1] = np.random.normal(0., muBm4ClDmo.TransSigma(), None)
        Rmu[1][2] = 0.

        Pmu[1][0] = np.random.normal(0., muBm4ClDmo.TransMmtm(), None)
        Pmu[1][1] = np.random.normal(0., muBm4ClDmo.TransMmtm(), None)
        delE      = np.random.normal(0., muBm4ClDmo.EnergySigma(), None)
        Pmu[0]    = np.sqrt(mc.mass()**2 + muBm4ClDmo.MmtmBar()**2) + delE
        Ptot2     = Pmu[0]**2 - mc.mass()**2
        Pmu[1][2] = np.sqrt(Ptot2 - Pmu[1][0]**2 - Pmu[1][1]**2)

        self._Rmu = Rmu
        self._Pmu = Pmu
        
    def __repr__(self):
        return "mu4ClDmo()"

    def __str__(self):
        return "mu4ClDmo: version=%g \n" \
               "    Rmu=[%g, (%g, %g, %g)] ns, mm \n" \
               "    Pmu=[%g, (%g, %g, %g)] MeV, MeV/c \n" \
               % \
               (self.CdVrsn(), \
                self.Rmu()[0], self.Rmu()[1][0],self.Rmu()[1][1],self.Rmu()[1][2], \
                self.Pmu()[0], self.Pmu()[1][0],self.Pmu()[1][1],self.Pmu()[1][2] \
                )

#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)

    def CdVrsn(self):
        return 1.0

    def Rmu(self):
        return self._Rmu

    def Pmu(self):
        return self._Pmu
