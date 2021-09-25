#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class PionEventInstance:
============================

  Generates single muon production event along the nuSTORM production
  straight.

  Class attributes:
  -----------------
  __MuonDecay: muon decay class
  --np       : numpy class
      
  Instance attributes:
  --------------------
  _ppi      : Pion momentum: i/p argument at instance creation
  _pmu      : Muon momentum; 
  _TrcSpcCrd: Trace space (s, x, y, z, x', y') in numpy array at 
              point of decay
  _pmuGen   : Generated pion momentum
  _P_mu     : Muon 4 momentum: (E, array(px, py, pz)), GeV
  _P_numu   : Muon-neuutrino 4 momentum: (E, array(px, py, pz)), GeV
  _tpi      : Time of pion decay
    
  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates decay instance
      __repr__ : One liner with call.
      __str__  : Dump of values of decay

  Get/set methods:
    getpppi           : Returns pion momentum as real
    getpmu            : Returns muon momentum as real
    getTraceSpaceCoord: Returns trace space: (s, x, y, z, x', y') (m)
    getppiGen         : Returns generated pion momentum (GeV)
    getmu4mmtm        : Returns muon 4 momentum: (E, array(px, py, pz)), GeV
    getnumu4mmtm      : Returns muon-neutrino 4 momentum: (E, array(px, py, pz)), GeV
    gettpi            : Returns time of pion decay
    
  General methods:
    CreateNeutrinos      : Manager for neutrino decay, returns z (m) of decay (real), P_e, P_nue, P_numu, RestFrame
                           Restframe contains a dump of the instance attributes of the MuonDecay class
    GenerateDcyPhaseSpace: Generates trace-space position of decay
    GenerateLongiPos     : Returns s, z of decay
    Boost2nuSTORM        : Boots to nuSTORM rest frame -- i.e. boost to ppi
    RotnBoost            : Operator; rotates and boosts rest-frame coordinates to nuSTORM frame

Created on Tue 30Mar21;02:26: Version history:
----------------------------------------------
 1.0: 30Mar21: First implementation - checked on 12 May that the decay length for a 6 GeV pion beam
 is as expected, by fitting the distribution of decay lengths

@author: kennethlong
@auhtor: PaulKyberd
"""

from copy import deepcopy
import nuSTORMPrdStrght as nuPrdStrt
import PionDecay as PionDecay
import MuonConst as MuonConst
import PionConst as PionConst
import numpy as np

muCnst = MuonConst.MuonConst()
piCnst = PionConst.PionConst()
nuStrt = nuPrdStrt.nuSTORMPrdStrght('11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
    
class PionEventInstance:

    __mumass = muCnst.mass()/1000.
    __pimass = piCnst.mass()/1000.
    __sol    = muCnst.SoL()

    __Debug  = False
    
#--------  "Built-in methods":
    def __init__(self, ppi):

        self._ppi = ppi
        self._phi=-10.0
        self._costheta=-10.0
        self._TrcSpcCrd, self._ppiGen, self._P_mu, self._P_numu, self._tpi = self.CreateMuon()

        return

    def __repr__(self):
        return "PionEventInstance(pmu)"

    def __str__(self):
        return "PionEventInstance: tpi (ns or s?) = %g, p_mu (GeV) = %g, s (m) = %g, z (m) = %g, generated momentum=%g, \r\n \
                P_mu (%g, [%g, %g, %g]), \r\n \
                P_numu (%g, [%g, %g, %g]), \r\n" % \
            (self._tpi, self._ppi, self._TrcSpcCrd[0], self._TrcSpcCrd[3], self._ppiGen, \
             self._P_mu[0], self._P_mu[1][0],self._P_mu[1][1],self._P_mu[1][2], \
             self._P_numu[0], self._P_numu[1][0],self._P_numu[1][1],self._P_numu[1][2] )
             
    
#--------  Generation of neutrino-creation event:
#.. Manager:
    def CreateMuon(self):
#        print ("Starting create muon")
        
        PrdStrghtLngth = 180.0
#        #.. Prepare--get muon decay instance in pion rest frame:
        z = 2.* PrdStrghtLngth
        Dcy = 0
        if PionEventInstance.__Debug:
            print("PionEventInstance.CreatePion: find valid decay")
        while z > PrdStrghtLngth:
            if isinstance(Dcy, PionDecay.PionDecay):
                del Dcy
            
            Ppi0 = self.getppi()
#            print ("Ppi0 ", Ppi0)
#   Comment out this line for distribution
#            Ppi = Ppi0
#    Comment out this line and we get single energy
            Ppi = nuStrt.GenerateMmtm(Ppi0)
#            print ("Ppi ", Ppi)
            Epi   = np.sqrt(Ppi**2 + PionEventInstance.__pimass**2)
            beta  = Ppi / Epi
            gamma = Epi / PionEventInstance.__pimass
            v    = beta * PionEventInstance.__sol
            Tmax = nuStrt.ProdStrghtLen() / (gamma * v)

            Dcy = PionDecay.PionDecay(Tmax=Tmax)
            tpi = Dcy.getLifetime()
            self._phi = Dcy.getphi()
            self._costheta = Dcy.getcostheta()
#            print("PionEventInstance.CreatePion: Dcy ", Dcy)
            DcyCoord, ppiGen = self.GenerateDcyPhaseSpace(Dcy,Ppi)
#            print("PionEventInstance.CreatePion: ppiGen ", ppiGen)
#            print("PionEventInstance.CreatePion: DcyCoord ", DcyCoord)
            z = DcyCoord[3]
        if z > PrdStrghtLngth:
            print("PionEventInstance.CreatePion Alarm:", z)
        if PionEventInstance.__Debug:
            print("PionEventInstance.CreatePion: decay at z =", z)
            print("----> Dcy:", Dcy.__str__())
            
        #.. Boost to nuSTORM frame:
        if PionEventInstance.__Debug:
            print("PionEventInstance.CreateMuon: rotate and boost to nuSTORM rest frame:")
        P_mu, P_numu = self.Boost2nuSTORM(Dcy,Ppi)

        del Dcy

#        print ("Ending create muon")
        return DcyCoord, ppiGen, P_mu, P_numu, tpi

#.. Trace space coordinate generation: array(s, x, y, z, x', y')
    def GenerateDcyPhaseSpace(self, Dcy, Ppi):
 #       print ("generate phase space start")
        coord = np.array([0., 0., 0., 0., 0., 0.])
        #.. longitudinal position, "s", z:
#        coord[0] = self.GenerateLongiPos(Dcy, Ppi)
# using  s = (beta*c) * (gamma*t) = p/E * c * E/m * lifetime = p*c*lifetime/m 
        coord[0] =Ppi*PionEventInstance.__sol*Dcy.getLifetime()/PionEventInstance.__pimass
        # [3] is the z position which is calculated from the decay distance and information about the
        # beam trajectory. At present set to the valus of s
        coord[3] = coord[0]
# here we set the Beam Twiss parameters and Energy spread 
        x, y, xp, yp = nuStrt.GenerateTrans(coord[0])
        coord[1] = x
        coord[2] = y
        coord[4] = xp
        coord[5] = yp

#        print ("generate phase space end")
        return coord, Ppi

#.. Boost from muon rest frame to nuSTORM frame:
    def Boost2nuSTORM(self, Dcy, PPi):
        ''' Present approximation is muon propagates along z axis, so, boost only 
            in preparation for later, include rotation matrix to transform from
            nustorm frame to frame with z axis along muon momentum and back '''

#        print ("PPi on Start boost to nuStorm", PPi)

       
        EPi = np.sqrt(PPi**2 + PionEventInstance.__pimass**2)
        beta   = PPi / EPi
        gamma  = EPi / PionEventInstance.__pimass

        R    = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        Rinv = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        if PionEventInstance.__Debug:
            print("PionEventInstance.Boost2nuSTORM: boost parameters:")
            print("----> PPi, EPi, beta, gamma:", PPi, EPi, beta, gamma)    # OK with P_mu=8

        # Treat decay components:
        P_mu    = Dcy.get4vmu()
        if PionEventInstance.__Debug:
            print("PionEventInstance.Boost2nuSTORM: rest frame P_mu:")
            print("----> Dcy.get4vmu:", Dcy.get4vmu())
            print("----> P_mu (MeV)", P_mu)
        P_mu[0] = P_mu[0]/1000.
        P_mu[1] = P_mu[1]/1000.
        if PionEventInstance.__Debug:
            print("----> P_mu (GeV):", P_mu)
            print("----> Dcy.get4vmu:", Dcy.get4vmu())
        P_mu    = self.RotnBoost(P_mu, R, Rinv, gamma, beta)
        if PionEventInstance.__Debug:
            print("PionEventInstance.Boost2nuSTORM: nuSTORM frame P_mu:")
            print("----> P_mu (GeV):", P_mu)
            print("----> Dcy.get4vmu:", Dcy.get4vmu())
        
        P_numu    = Dcy.get4vnumu()
        P_numu[0] = P_numu[0]/1000.
        P_numu[1] = P_numu[1]/1000.
        P_numu    = self.RotnBoost(P_numu, R, Rinv, gamma, beta)
        
 #       print ("End boost to nuStorm")

        return P_mu, P_numu
    
    def RotnBoost(self, P, R, Rinv, gamma, beta):
        if PionEventInstance.__Debug:
            print("PionEventInstance.RotnBoos:")
            print("----> P, R, Rinv, gamma, beta:", P, "\n", \
                  R, "\n", Rinv, "\n", gamma, beta)
            
  #      print ("P[1] is ", P[1])
  #      print ("P is ", P)
        p3 = np.dot(R, P[1])
 #       print("    p3 is ", p3)
        
        Ec = P[0]
        Pc = p3[2]
#        print("    Ec is ", Ec, "   Pc is ", Pc)
        
        Ef = gamma * (Ec + beta * Pc)
        Pf = gamma * (Pc + beta * Ec)
#        print("    Ef is ", Ef, "   Pf is ", Pf)
        
        p3[2] = Pf
        p3    = np.dot(Rinv, p3)
#        print("    p3 is ", p3)

        Po    = [0., np.array([0., 0., 0.])]
        Po[0] = Ef
        Po[1] = p3
#        print("    Po is ", Po)
        
        return Po
    
#--------  get/set methods:
    def gettpi(self):
        return deepcopy(self._tpi)
        
    def getppi(self):
        return deepcopy(self._ppi)

    def getphi(self):
        return deepcopy(self._phi)

    def getcostheta(self):
        return deepcopy(self._costheta)

    def getTraceSpaceCoord(self):
        return deepcopy(self._TrcSpcCrd)

    def getmu4mmtm(self):
        return deepcopy(self._P_mu)

    def getnumu4mmtm(self):
        return deepcopy(self._P_numu)

    def getppiGen(self):
        return deepcopy(self._ppiGen)
