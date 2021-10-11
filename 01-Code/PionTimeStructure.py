#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class PionTimeStr:
==================

  Generates time distributions of pions at the decay straight entrance.

#NEED TO ADJUST!!!!
  Dependencies:
   - numpy, Simulation, math, PionConst

#NEED TO ADJUST!!!!
  Class attributes:
  -----------------
  __piCnst: Instance of PionConst class.

#NEED TO ADJUST!!!!
  Instance attributes:
  --------------------
  _Lifetime : Time to decay (s) multiplied by the speed of light.  This is
              is for "this particular pion"
  _v_numu   : Muon-neutrino 4-vector (E, array(px, py, pz)); MeV
  _v_mu        : Muon 4-vector (E, array(px, py, pz)); MeV
  _costheta : Cosine of angle between the pion and the decay muon in the pion rest frame
  _phi      : Angle polar angle of the decay between the pion and the decay muon in the pion rest
              frame - measured with respect to the x axis of the nuStorm co-ordinate system

  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates decay instance
      __repr__ : One liner with call.
      __str__  : Dump of values of decay

  Get/set methods:
    getCrossingtime: Returns time (s) of this pion instance crossing the entrance to production straight


#NEED TO ADJUST!!!!
  Pion-decay methods:
    GenerateLifetime: Generates lifetime of this instance.  Returns
                      lifetime (float). Units s
    decaypion         : Generates a particular muon decay; calls each of
                      following methods in turn.  Returns 32 4-vectors,
                      v_mu, v_numu.  v_i = [Energy, array(px, py, px)]
                      (units MeV), and two floats, costheta, phi
    get3vectors     : Generates 3-vector momenta (MeV).  muon direction taken
                      as positive z direction.  Returns two float array objects,
                      p_mu, p_numu; p_i = array(px, py, pz), and two floats
                      costheta and phi
    ranCoor         : Applies random 3D rotation and returns 4-vectors:
                      v_mu, v_numu.  v_i = [Energy, array(px, py, px)]
                      Units MeV.

Created on 17Jul21; Version history:
----------------------------------------------
 1.0: 17Jul21: First implementation - based on standalone code and PionDecay

@author: MarvinPfaff
"""

import numpy as np
import PionConst as pC

piCnst = pC.PionConst()


class PionTimeStructure:

    __pimass = piCnst.mass()/1000.
    __sol = piCnst.SoL()

    __ltransfer = 50.

    __Debug = False

#--------  "Built-in methods":
    def __init__(self, **kwargs):

        ppi = kwargs.get('ppi', 8.00)           #pion momentum
        t0 = kwargs.get('t0', 0.)               #time inside bunch when first pion leaves target

        self._Crossingtime = self.GenerateCrossingtime(ppi = ppi, t0=t0)

        return

    def __repr__(self):
        return "PionTimeStructure()"

    def __str__(self):
        return "PionTimeStructure: Crossingtime=%g,      \
               " % (self._Crossingtime)

#--------  "Dynamic methods"; individual lifetime, energies, and angles
    def GenerateCrossingtime(self, **kwargs):
        Epi  = np.sqrt(ppi**2 + PionTimeStructure.__pimass**2)
        beta = ppi / Epi
        v    = beta * PionTimeStructure.__sol   #depending on implementation, v can also be handed over directly

        t = PionTimeStructure.__ltransfer / v + t0

        return t


#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)
    def getCrossingtime(self):
        return self._Crossingtime
