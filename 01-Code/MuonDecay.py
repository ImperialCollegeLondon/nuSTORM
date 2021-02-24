#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class MuonDecay:
================

  Generates decay distributions for muon decay at rest.  Kinematic 
  distributions are for muon decay at rest.

  Citations: 
   - G4 manual

  Dependencies:
   - numpy, Simulation, math, MuonConst

  Class attributes:
  -----------------
  __muCnst: Instance of MuonConst class.
      
  Instance attributes:
  --------------------
  _Lifetimre: Time to decay (s) multiplied by the speed of light.  This is
              is for "this particular muon"
  _v_e      : Electron 4-vector (E, array(px, py, pz)); MeV
  _v_nue    : Electron-neutrino 4-vector (E, array(px, py, pz)); MeV
  _v_numu   : Muon-neutrino 4-vector (E, array(px, py, pz)); MeV
  _costheta : Cosine of angle between electron and nue in muon rest frame
  _cosphi   : Cosine of angle between electron and numu in muon rest frame

  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates decay instance
      __repr__ : One liner with call.
      __str__  : Dump of values of decay

  Get/set methods:
    getLifetime: Returns lifetime (s) of this decay instance
    get4ve     : Returns 4-vector of electron (MeV)
    get4vnue   : Returns 4-vector of electron neutrino (MeV)
    get4vnumu  : Returns 4-vector of muon neutrino (MeV)
    getcostheta: Returns cosine of angle between electron and 
                 electron-neutrino in muon rest frame
    getcosphi  : Returns cosine of angle between electron and 
                 muon-neutrino in muon rest frame
  
  Muon-decay methods:
    GenerateLifetime: Generates lifetime of this instance.  Returns
                      lifetime (float). Units s
    decaymuon       : Generates a particular muon decay; calls each of 
                      following methods in turn.  Returns 3 4-vectors,
                      v_e, v_nue, v_numu.  v_i = [Energy, array(px, py, px)] 
                      (units MeV), and two floats, costheta, cosphi
    GenerateScldE   : Generates scaled energies of electron, nu_e, and nu_mu.  
                      Returns three floats, (f_e, f_nue, f_numu).  Units m_mu/2.
    get3vectors     : Generates 3-vector momenta (MeV).  Electron direction taken
                      as positive z direction.  Returns three floar array objects,
                      p_e, p_mue, p_numu; p_i = array(px, py, pz), and two floats
                      costheta and cosphi
    ranCoor         : Applies random 3D rotation and returns 4-vectors:
                      v_e, v_nue, v_numu.  v_i = [Energy, array(px, py, px)]
                      Units MeV.

Created on Sat 09Jan21;17:18: Version history:
----------------------------------------------
 1.0: 09Jan21: First implementation

@author: kennethlong
"""

from copy import deepcopy 
#from Simulation import *
#from Simulation import getRandom
import Simulation as Simu
import numpy as np
import math as mth
import MuonConst as muConst

muCnst = muConst.MuonConst()

class MuonDecay:
    
    __Debug = False
    
#--------  "Built-in methods":
    def __init__(self):

        self._Lifetime = self.GenerateLifetime()

        v_e, v_nue, v_numu, costheta, cosphi = self.decaymuon()

        self._v_e = v_e
        self._v_nue = v_nue
        self._v_numu = v_numu
        
        self._costheta = costheta
        self._cosphi   = cosphi
        
        return

    def __repr__(self):
        return "MuonDecay()"

    def __str__(self):
        return "MuonDecay: Lifetime=%g, \r\n \
                           v_e=(%g, [%g, %g, %g]), \r\n \
                         v_nue=(%g, [%g, %g, %g]), \r\n \
                        v_numu=(%g, [%g, %g, %g]),      \
               " % (self._Lifetime,                                                       \
                     self._v_e[0], self._v_e[1][0], self._v_e[1][1], self._v_e[1][2],              \
                     self._v_nue[0], self._v_nue[1][0], self._v_nue[1][1], self._v_nue[1][2],      \
                     self._v_numu[0], self._v_numu[1][0], self._v_numu[1][1], self._v_numu[1][2])

#--------  "Dynamic methods"; individual lifetime, energies, and angles
    def GenerateLifetime(self):
        ran = Simu.getRandom()
        lt = -mth.log(1.-ran) * muCnst.lifetime()
        return lt

    def GenerateScldE(self):
#----  See P. 3, my "Notes and calcs"
        f_nue = 0.5
        f_e   = 0.
        while f_nue < (1. - f_e):
#.. fractional electron energy:
            f_e = Simu.getRandom()
            
#.. fractional electron-neutrino energy:
            gam = 1.
            x   = 1.
            while gam > x*(1.-x):
                x   = Simu.getRandom()
                gam = Simu.getRandom()
            f_nue = x

        f_numu = 2. - f_e - f_nue
        
        return f_e, f_nue, f_numu

    def get3vectors(self, f_e, f_nue, f_numu):
#----  See P. 3, my "Notes and calcs"
        #print("f_e, f_nue, f_numu:", f_e, f_nue, f_numu)
#.. electron neutrino
        costheta = 1. - 2.*( 1./f_e + 1/f_nue - 1/(f_e*f_nue) )
        sintheta = mth.sqrt(1. - costheta**2)
        #print("theta:", sintheta, costheta)
        
#.. muon neutrino
        cosphi = -(f_e + f_nue*costheta) / f_numu
        sinphi = -f_nue*sintheta / f_numu
        diff = abs(1. - (cosphi**2 + sinphi**2))

        theta = mth.acos(costheta) * 180./mth.pi
        phi   = mth.acos(cosphi)   * 180./mth.pi
        #print("phi:", sinphi, cosphi)

#.. Three vectors:
        p_e    = np.array([             0., 0., f_e            ])
        p_nue  = np.array([ f_nue*sintheta, 0., f_nue*costheta ])
        p_numu = np.array([  f_numu*sinphi, 0., f_numu*cosphi  ])
        
        return p_e, p_nue, p_numu, costheta, cosphi

    def ranCoor(self, p_e, p_nue, p_numu):
#.. Rotation angles
        alpha  = Simu.getRandom() * 2.*mth.pi
        calpha = mth.cos(alpha)
        salpha = mth.sin(alpha)

        cbeta = -1. + 2.*Simu.getRandom()
        sbeta = mth.sqrt(1. - cbeta**2)

        gamma = Simu.getRandom() * 2.*mth.pi
        cgamma = mth.cos(gamma)
        sgamma = mth.sin(gamma)

        #print("--------  --------  --------  --------  -------- ")
        #print("p_e", p_e)
        #print("p_nue", p_nue)
        #print("p_numu", p_numu)

#.. Rotation angles
        Ra = np.array([          \
             [calpha   , -salpha,      0.      ], \
             [salpha   ,  calpha,      0.      ], \
             [0.       , 0.,      1.      ] \
                                   ])
        #print("Ra:", Ra)
        Rb = np.array([          \
             [cbeta    , 0.,      -sbeta     ], \
             [0.       , 1.,      0.      ], \
             [sbeta    , 0.,      cbeta      ] \
                                   ])
        #print("Rb:", Rb)
        Rc = np.array([          \
             [cgamma   , -sgamma,      0.      ], \
             [sgamma   ,  cgamma,      0.      ], \
             [0.       , 0.,      1.      ] \
                                   ])
        #print("Rc:", Rc)

        Rr = np.dot(Ra, Rb)
        #print("Rr:", Rr)
        Rr = np.dot(Rr, Rc)
        #print("Rr:", Rr)
        
#.. Do rotation:
        p_e1    = np.dot(Rr, p_e)
        p_nue1  = np.dot(Rr, p_nue)
        p_numu1 = np.dot(Rr, p_numu)
        #print("p_e1", p_e1)
        #print("p_nue1", p_nue1)
        #print("p_numu1", p_numu1)

        return p_e1, p_nue1, p_numu1

    def decaymuon(self):
#.. Get scaled muon eneries:        
        f_e, f_nue, f_numu = self.GenerateScldE()
        
#.. Get scaled 3-vectors:
        s_e, s_nue, s_numu, costheta, cosphi = self.get3vectors(f_e, f_nue, f_numu)

#.. Rotate to arbitrary axis orientation:
        p_e, p_nue, p_numu = self.ranCoor(s_e, s_nue, s_numu)

#.. Scale to GeV and make 4-vectors:
        mo2    = muCnst.mass() / 2.
        
        f_e    = f_e    * mo2
        p_e    = p_e    * mo2
        
        f_nue  = f_nue  * mo2
        p_nue  = p_nue  * mo2
        
        f_numu = f_numu * mo2
        p_numu = p_numu * mo2

        v_e    = [f_e,    p_e]
        v_nue  = [f_nue,  p_nue]
        v_numu = [f_numu, p_numu]

        return v_e, v_nue, v_numu, costheta, cosphi

#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)
    def getLifetime(self):
        return self._Lifetime

    def get4ve(self):
        return deepcopy(self._v_e)

    def get4vnue(self):
        return deepcopy(self._v_nue)

    def get4vnumu(self):
        return deepcopy(self._v_numu)

    def getcostheta(self):
        return deepcopy(self._costheta)

    def getcosphi(self):
        return deepcopy(self._cosphi)
    
