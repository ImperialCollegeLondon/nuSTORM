#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class PionDecay:
================

  Generates decay distributions for pion decay at rest.  Kinematic 
  distributions are for pion decay at rest.

  Dependencies:
   - numpy, Simulation, math, PionConst

  Class attributes:
  -----------------
  __piCnst: Instance of PionConst class.
      
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
    getLifetime: Returns lifetime (s) of this decay instance
    get4vnumu  : Returns 4-vector of muon neutrino (MeV)
    get4mu     : Returns 4-vector of muon (MeV)
    getcostheta: Returns cosine of angle between muon and the parent pion in the pion rest frame
    getphi     : The polar  angle of the muon in the pion rest frame, wrt the x axis of nuStorm

  
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

Created on 26March21; Version history:
----------------------------------------------
 1.0: 26Mar21: First implementation - based on MuonDecay

@author: kennethlong and PaulKyberd
"""

from copy import deepcopy 
import Simulation as Simu
import numpy as np
import math as mth
import PionConst as pC
import MuonConst as mC

piCnst = pC.PionConst()
muCnst = mC.MuonConst()

class PionDecay:

    __Debug = False    
    
#--------  "Built-in methods":
    def __init__(self, **kwargs):

        Tmax = kwargs.get('Tmax', float('inf'))
        ppi = kwargs.get('ppi', 8.00)
        
        self._Lifetime = self.GenerateLifetime(Tmax=Tmax)

        v_mu, v_numu, costheta, phi = self.decaypion()

        self._v_mu = v_mu
        self._v_numu = v_numu
        
        self._costheta = costheta
        self._phi   = phi
        
        return

    def __repr__(self):
        return "PionDecay()"

    def __str__(self):
        return "PionDecay: Lifetime=%g, \r\n \
                           v_mu=(%g, [%g, %g, %g]), \r\n \
                           v_numu=(%g, [%g, %g, %g]),      \
               " % (self._Lifetime,                                                       \
                     self._v_mu[0], self._v_mu[1][0], self._v_mu[1][1], self._v_mu[1][2],              \
                     self._v_numu[0], self._v_numu[1][0], self._v_numu[1][1], self._v_numu[1][2])

#--------  "Dynamic methods"; individual lifetime, energies, and angles
    def GenerateLifetime(self, **kwargs):
        Tmax = kwargs.get('Tmax', float('inf'))
        Gmx = 1. - mth.exp( -Tmax / piCnst.lifetime() )
        ran = Simu.getRandom() * Gmx
        lt = -mth.log(1.-ran) * piCnst.lifetime()
        return lt


#--------  Random rotation flat in phi and cos Theta
    def ranCoor(self, p_mu, p_numu):
#.. Rotation angles
        phi = Simu.getRandom() * 2.*mth.pi
        cPhi = mth.cos(phi)
        sPhi = mth.sin(phi)

        cTheta = -1. + 2.*Simu.getRandom()
        sTheta = mth.sqrt(1. - cTheta**2)
        theta = mth.acos(cTheta)

#.. Rotation angles
        Ra = np.array([          \
             [cPhi   , -sPhi,    0.      ], \
             [sPhi   ,  cPhi,    0.      ], \
             [0.     , 0.,       1.      ] \
                                ])
#        print("Ra:", Ra)
        Rb = np.array([          \
             [cTheta    , 0.,    -sTheta     ], \
             [0.        , 1.,    0.          ], \
             [sTheta    , 0.,    cTheta      ] \
                                   ])
#        print("Rb:", Rb)

        Rr = np.dot(Ra, Rb)
#        print("Rr:", Rr)
        
#.. Do rotation:
        p_mu1  = np.dot(Rr, p_mu)
        p_numu1 = np.dot(Rr, p_numu)

        return  p_mu1, p_numu1, cTheta, phi

    def decaypion(self):
#.. Calculate 3-vectors - muon is going forward:
        f_mu = [0.0, 0.0, 29.79]
        f_numu = [0.0, 0.0, -29.79]
#        print ("muon 3 vectors ", f_mu, "    numu 3 vector ", f_numu)

#.. Rotate to arbitrary axis orientation:
#        print ("rotating vectors")
        
        p_mu, p_numu, cosTheta, phi = self.ranCoor(f_mu, f_numu)

#        print ("back from rotation ... p_mu ", p_mu, "    p_numu  ", p_numu)

        E_mu = mth.sqrt(muCnst.mass()*muCnst.mass() + p_mu[0]*p_mu[0] + p_mu[1]*p_mu[1] + p_mu[2]*p_mu[2])
        E_numu = 29.79
#  Energy conservation check - put the check in the debug check when happy
        if self.__Debug == True:
            ETot = E_mu + E_numu
            print ("Check energy conservation Emu = ", E_mu, "    Enumu = ", E_numu, "    ETotal = ", ETot)
        ETot = E_mu + E_numu
        if abs(ETot-piCnst.mass()) > 0.01:
            print ("Energy not conserved, delta E is ", ETot-piCnst.mass())
        
        v_mu   = [E_mu,  p_mu]
        v_numu = [E_numu, p_numu]

        return v_mu, v_numu, cosTheta, phi

#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)
    def getLifetime(self):
        return self._Lifetime

    def get4vmu(self):
        return deepcopy(self._v_mu)

    def get4vnumu(self):
        return deepcopy(self._v_numu)

    def getcostheta(self):
        return deepcopy(self._costheta)

    def getphi(self):
        return deepcopy(self._phi)
    
