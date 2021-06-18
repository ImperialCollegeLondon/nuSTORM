#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class NeutrinoEventInstance:
============================

  Generates single neutrino-production event along the nuSTORM production
  straight.

  Class attributes:
  -----------------
  __MuonDecay: muon decay class
  --np       : numpy class
      
  Instance attributes:
  --------------------
  _pmu      : Muon momentum; i/p argument at instance creation
  _ct       : Decay time * speed-of-light
  _TrcSpcCrd: Trace space (s, x, y, z, x', y') in numpy array at 
              point of decay
  _pmuGen   : Generated muon momentum
  _P_e      : Electron 4 momentum: (E, array(px, py, pz)), GeV
  _P_nue    : Electron-neutrino 4 momentum: (E, array(px, py, pz)), GeV
  _P_numu   : Muon-neuutrino 4 momentum: (E, array(px, py, pz)), GeV
    
  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates decay instance
      __repr__ : One liner with call.
      __str__  : Dump of values of decay

  Get/set methods:
    getpmu            : Returns muon momentum as real
    getTraceSpaceCoord: Returns trace space: (s, x, y, z, x', y') (m)
    getpmuGen         : Returns generated muon momentum (GeV)
    gete4mmtm         : Returns electron 4 momentum: (E, array(px, py, pz)), GeV
    getnue4mmtm       : Returns electron-neutrino 4 momentum: (E, array(px, py, pz)), GeV
    getnumu4mmtm      : Returns muon-neutrino 4 momentum: (E, array(px, py, pz)), GeV
  
  General methods:
    CreateNeutrinos      : Manager for neutrino decay, returns z (m) of decay (real), P_e, P_nue, P_numu, RestFrame
                           Restframe contains a dump of the instance attributes of the MuonDecay class
    GenerateDcyPhaseSpace: Generates trace-space position of decay
    GenerateLongiPos     : Returns s, z of decay
    Boost2nuSTORM        : Boots to nuSTORM rest frame -- i.e. boost to pmu
    RotnBoost            : Operator; rotates and boosts rest-frame coordinates to nuSTORM frame

Created on Sat 16Jan21;02:26: Version history:
----------------------------------------------
 1.2: 18Jun21: Pass nuSTORM production straight parameter file name
 1.1: 03Apr21: Fix error in relativistic treatment muon lifetime
 1.0: 16Jan21: First implementation

@author: kennethlong
"""

from copy import deepcopy
import nuSTORMPrdStrght as nuPrdStrt
import MuonDecay as MuonDecay
import MuonConst as MuonConst
import numpy as np

muCnst = MuonConst.MuonConst()
    
class NeutrinoEventInstance:

    __mumass = muCnst.mass()/1000.
    __sol    = muCnst.SoL()

    __Debug  = False

#--------  "Built-in methods":
    def __init__(self, pmu=5., filename=None):

        nuStrt = nuPrdStrt.nuSTORMPrdStrght(filename)

        self._pmu = pmu
        self._ct, self._TrcSpcCrd, self._pmuGen, self._P_e, self._P_nue, self._P_numu = self.CreateNeutrinos(nuStrt)

        return

    def __repr__(self):
        return "NeutrinoEventInstance(pmu)"

    def __str__(self):
        return "NeutrinoEventInstance: p_mu (GeV) = %g, s (m) = %g, z (m) = %g, generated momentum=%g, \r\n \
                P_e (%g, [%g,, %g, %g]), \r\n \
                P_nue (%g, [%g,, %g, %g]), \r\n \
                P_numu (%g, [%g,, %g, %g]), \r\n" % \
            (self._pmu, self._TrcSpcCrd[0], self._TrcSpcCrd[3], self._pmuGen, \
             self._P_e[0], self._P_e[1][0], self._P_e[1][1], self._P_e[1][2], \
             self._P_nue[0], self._P_nue[1][0],self._P_nue[1][1],self._P_nue[1][2], \
             self._P_numu[0], self._P_numu[1][0],self._P_numu[1][1],self._P_numu[1][2] )
    
#--------  Generation of neutrino-creation event:
#.. Manager:
    def CreateNeutrinos(self, nuStrt):
        PrdStrghtLngth = nuStrt.ProdStrghtLen()
        #.. Prepare--get neutrino decay instance in muon rest frame:
        z = 2.* PrdStrghtLngth
        Dcy = 0
        if NeutrinoEventInstance.__Debug:
            print("NeutrinoEventInstance.CreateNeutrinos: find valid decay")
        while z > PrdStrghtLngth:
            if isinstance(Dcy, MuonDecay.MuonDecay):
                del Dcy
            
            Pmu0 = self.getpmu()
            Pmu  = nuStrt.GenerateMmtm(Pmu0)
            Emu   = np.sqrt(Pmu**2 + NeutrinoEventInstance.__mumass**2)
            beta  = Pmu / Emu
            gamma = Emu / NeutrinoEventInstance.__mumass
            v    = beta * NeutrinoEventInstance.__sol
            Tmax = nuStrt.ProdStrghtLen() / (gamma * v)

            Dcy = MuonDecay.MuonDecay(Tmax=Tmax)
            DcyCoord, pmuGen = self.GenerateDcyPhaseSpace(Dcy, Pmu, nuStrt)
            z  = DcyCoord[3]
            ct = Dcy.getLifetime()

        if z > PrdStrghtLngth:
            print("NeutrinoEvenInstance.CreateNeutrinos Alarm:", z)
        if NeutrinoEventInstance.__Debug:
            print("NeutrinoEventInstance.CreateNeutrinos: decay at z =", z)
            print("----> Dcy:", Dcy.__str__())
            
        #.. Boost to nuSTORM frame:
        if NeutrinoEventInstance.__Debug:
            print("NeutrinoEventInstance.CreateNeutrinos: rotate and boost to nuSTORM rest frame:")
        P_e, P_nue, P_numu = self.Boost2nuSTORM(Dcy)
        if NeutrinoEventInstance.__Debug:
            print("----> P_e   :", P_e)
            print("----> P_nue :", P_nue)
            print("----> P_numu:", P_numu)

        del Dcy
        return ct, DcyCoord, pmuGen, P_e, P_nue, P_numu

#.. Trace space coordinate generation: array(s, x, y, z, x', y')
    def GenerateDcyPhaseSpace(self, Dcy, Pmu, nuStrt):
        coord = np.array([0., 0., 0., 0., 0., 0.])

        #.. longitudinal position, "s", z:
        coord[0] = self.GenerateLongiPos(Dcy, Pmu)
        coord[3] = nuStrt.Calculatez(coord[0])

        x, y, xp, yp = nuStrt.GenerateTrans(coord[0])
        coord[1] = x
        coord[2] = y
        coord[4] = xp
        coord[5] = yp

        return coord, Pmu

#.. Trace space coordinate generation:
    def GenerateLongiPos(self, Dcy, Pmu):
        Emu   = np.sqrt(Pmu**2 + NeutrinoEventInstance.__mumass**2)
        beta  = Pmu / Emu
        gamma = Emu / NeutrinoEventInstance.__mumass
        
        v    = beta * NeutrinoEventInstance.__sol

        s   = v * gamma * Dcy.getLifetime()
        
        return s

#.. Boost from muon rest frame to nuSTORM frame:
    def Boost2nuSTORM(self, Dcy):
        ''' Present approximation is muon propagates along z axis, so, boost only 
            in preparation for later, include rotation matrix to transform from
            nustorm frame to frame with z axis along muon momentum and back '''

        Pmu = self.getpmu()
        Emu = np.sqrt(Pmu**2 + \
                                              NeutrinoEventInstance.__mumass**2)
        beta   = Pmu / Emu
        gamma  = Emu / NeutrinoEventInstance.__mumass

        R    = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        Rinv = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        if NeutrinoEventInstance.__Debug:
            print("NeutrinoEventInstance.Boost2nuSTORM: boost parameters:")
            print("----> Pmu, Emu, beta, gamma:", Pmu, Emu, beta, gamma)

        # Treat decay components:
        P_e    = Dcy.get4ve()
        if NeutrinoEventInstance.__Debug:
            print("NeutrinoEventInstance.Boost2nuSTORM: rest frame P_e:")
            print("----> Dcy.get4ve:", Dcy.get4ve())
            print("----> P_e (MeV)", P_e)
        P_e[0] = P_e[0]/1000.
        P_e[1] = P_e[1]/1000.
        if NeutrinoEventInstance.__Debug:
            print("----> P_e (GeV):", P_e)
            print("----> Dcy.get4ve:", Dcy.get4ve())
        P_e    = self.RotnBoost(P_e, R, Rinv, gamma, beta)
        if NeutrinoEventInstance.__Debug:
            print("NeutrinoEventInstance.Boost2nuSTORM: nuSTORM frame P_e:")
            print("----> P_e (GeV):", P_e)
            print("----> Dcy.get4ve:", Dcy.get4ve())
        
        P_nue    = Dcy.get4vnue()
        P_nue[0] = P_nue[0]/1000.
        P_nue[1] = P_nue[1]/1000.
        P_nue    = self.RotnBoost(P_nue, R, Rinv, gamma, beta)
        
        P_numu    = Dcy.get4vnumu()
        P_numu[0] = P_numu[0]/1000.
        P_numu[1] = P_numu[1]/1000.
        P_numu    = self.RotnBoost(P_numu, R, Rinv, gamma, beta)
        
        return P_e, P_nue, P_numu
    
    def RotnBoost(self, P, R, Rinv, gamma, beta):
        if NeutrinoEventInstance.__Debug:
            print("NeutrinoEventInstance.RotnBoos:")
            print("----> P, R, Rinv, gamma, beta:", P, "\n", \
                  R, "\n", Rinv, "\n", gamma, beta)
            
        p3 = np.dot(R, P[1])
        
        Ec = P[0]
        Pc = p3[2]
        
        Ef = gamma * (Ec + beta * Pc)
        Pf = gamma * (Pc + beta * Ec)
        
        p3[2] = Pf
        p3    = np.dot(Rinv, p3)

        Po    = [0., np.array([0., 0., 0.])]
        Po[0] = Ef
        Po[1] = p3
        
        return Po
    
#--------  get/set methods:
    def getpmu(self):
        return deepcopy(self._pmu)

    def getTraceSpaceCoord(self):
        return deepcopy(self._TrcSpcCrd)

    def gete4mmtm(self):
        return deepcopy(self._P_e)

    def getnue4mmtm(self):
        return deepcopy(self._P_nue)

    def getnumu4mmtm(self):
        return deepcopy(self._P_numu)

    def getpmuGen(self):
        return deepcopy(self._pmuGen)
