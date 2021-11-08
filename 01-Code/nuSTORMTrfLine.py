#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class nuSTORMTrfLine:
=======================

  Defines the parameters of the nuSTORM transfer line and
  contains methods for the simulation of the transfer line optics.

  First implementation contains first order simulation in which
  s=z and representative values of emittance and beta are simulated
  for alpha=0 (see nuSIM-2021-01).

  Class attributes:
  -----------------
  __instance : Set on creation of first (and only) instance.

  Instance attributes:
  --------------------
  _filename        = Filename, including path, to cls file containing
                     nuSTORM transfer line specification
  _TrfLineParams   = Pandas dataframe containing specification in
                     _filename
  _TrfLineLen      = Production straight length (m)
  _piAcc           = Momentum acceptance (%)
  _epsilon         = Representative emittance/acceptance (pi mm rad)
  _beta            = Representative beta function (mm)
  _delT            = Bunch length at target (ns)

  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __new__ : Creates singleton class and prints version and values
                of transfer line parameters used.
      __repr__: One liner with call.
      __str__ : Dump of constants

  Get/set methods:
      CdVrsn()     : Returns code version number.
      TrfLineLen   : Get transfer line length (m)
      piAcc        : Get momentum acceptance (%)
      epsilon      : Get acceptance/emittance (pi mm rad)
      beta         : Get beta function (mm)
      delT         : Get bunch length (ns)


  Generation, calcution methods and utilities:
      GeneratePiMmtm : Generates momentum centred around p0 (input float).
                     Parabolic distribution generated.
      GenerateTime : Generates time to enter transfer line; first particle leaves at t=0.
                     Uniform distribution generated.
      GenerateTrans: Generate transverse phase space (x, y, xp, yp) given
                     representative emittance and beta.  Parabolic distributions
                     generated.
      GenerateDcyTime: Generate time in the nuSTORM frame after which pion decays.
      Calculatet   : Calculate t of particle at certain s_final, given initial s,
                     t_i and v.
      Calculatez   : Calculate z given s, p and transfer line
                     and production-straight parameters.
                     **This version: z=s** (m)
      Calculatepz  : Calculate pz from px and py.
      CalculateE   : Calculate E given p.
      Calculatev   : Calculate v given p.
      CalculateDcyPt : Calculate s of decay given gamma and v. It is also checked,
                     whether pion decays already in transfer line, given original s,
                     final s_final and intial t_i. If pion doesn't decay in transfer
                     line, s_dcy = -100. [m] is returned.
      GeneratePion : Generate two pion particle instances given an eventHistory instance,
                     runNum, eventNum, weight, plus initialization specific arguments.
                     These will be added to the eventHistory, which will then be returned.
                     The pion instance can be initialized in two different ways:
                     1. Pion is generated and initialized from scratch according to some
                        predefined distributions (so far - t: uniform, p: parabolic).
                        Additional argument needed is pion central momentum p0.
                     2. Pion is generated and initialized with input parameters. This
                        allows initialization with parameters from other simulations.
                        Additional arguments needed are x, y, px, py, pz and t.
                     The two pion instances are generated at target (s=-50m) and at
                     the entrance of the production straight (s=0m). It is also checked
                     whether the pion already decays in the transfer line. If so, the
                     pion at the production straight entrance has runNum=-1 and there
                     will be an additional pion, muon and neutrino at the decay point
                     added to the eventHistory.
      GenerateMuon : Generate a muon particle instance given runNum, eventNum, the
                     decay coordinates (DcyCoord: s,x,y,z), the muon momentum
                     (P_mu: px,py,pz), t and weight.
      GenerateNu   : Generate a neutrino particle instance given runNum, eventNum,
                     the decay coordinates (DcyCoord: s,x,y,z), the neutrino momentum
                     (P_numu: px,py,pz), t and weight.
      Boost2nuSTORM: Boots to nuSTORM rest frame -- i.e. boost to ppi.
                     This differs from PionEventInstance's methods in not taking an instance
                     of the PionDecay class as input but P_mu and P_numu in addition to ppi.
      RotnBoost    : Operator; rotates and boosts rest-frame coordinates to nuSTORM frame
      BoostBack2RestFrame: Boots to rest frame of pion decay -- i.e. inverting earlier boost
                     to nuSTORM rest frame.
      NegRotnBoost : Operator; rotates and boosts nuSTORM frame coordinates to rest-frame of pion
      getTrfLineParams: Read csv file and generate Pandas dataframe
      printParams  : Print all paramters


Created on Mo 01Nov21. Version history:
----------------------------------------
 1.0: 01Nov21: First implementation

@author: MarvinPfaff
based on nuSTORMPrdStrght class
"""

import pandas as pnds
import numpy as np
from copy import deepcopy
import PionConst as PionConst
import MuonConst as MuonConst
import PionDecay as PionDecay
import PionEventInstance as PionEventInstance
import particle as particle
import Simulation as Simu
import eventHistory as eventHistory

piCnst = PionConst.PionConst()
muCnst = MuonConst.MuonConst()


class nuSTORMTrfLine(object):
    __instance = None
    __pimass = piCnst.mass()/1000.
    __mumass = muCnst.mass()/1000.
    __piPDG = piCnst.pdgCode()
    __muPDG = muCnst.pdgCode()
    __piLifetime = piCnst.lifetime()
    __sol    = muCnst.SoL()

#--------  "Built-in methods":
    def __new__(cls, filename):
        if cls.__instance is None:
            print('nuSTORMTrfLine.__new__: creating the nuSTORMTrfLine object')
            cls.__instance = super(nuSTORMTrfLine, cls).__new__(cls)
            print('instance okay')
            cls._filename        = filename
            cls._TrfLineParams   = cls.GetTrfLineParams(filename)
            cls._TrfLineLen      = cls._TrfLineParams.iat[0,2]
            cls._piAcc            = cls._TrfLineParams.iat[1,2] / 100.
            cls._epsilon         = cls._TrfLineParams.iat[2,2]
            cls._beta            = cls._TrfLineParams.iat[3,2]
            cls._delT            = cls._TrfLineParams.iat[4,2]
            print('creation of new object: okay')
        return cls.__instance

    def __repr__(self):
        return "nuSTORMTrfLine()"

    def __str__(self):
        return "nuSTORMTrfLine: version=%g, length of transfer line=%g m, " \
               "momentum acceptance=%g%%, transverse acceptance=%g pi mm rad, beta=%g mm, bunch length=%g ns. \n" % \
               (self.CdVrsn(), self.TrfLineLen(), self.piAcc(), self.epsilon(), self.beta(), self.delT() )

#--------  Simulation methods:
    def GeneratePiMmtm(self,p0):
        p = -99.
        dp = p0 * self._piAcc
        p  = p0 + Simu.getParabolic(dp)
        return p

    def GenerateTime(self):
        t = Simu.getRandom() * self._delT * 10**(-9)
        return t

    def GenerateTrans(self,s):
        r  = np.sqrt(self._epsilon*self._beta) / 1000.
        rp = np.sqrt(self._epsilon/self._beta)
        x  = Simu.getParabolic(r)
        y  = Simu.getParabolic(r)
        xp = Simu.getParabolic(rp)
        yp = Simu.getParabolic(rp)
        return x, y, xp, yp

    def GenerateDcyTime(self,gamma):
        t = np.random.exponential(scale=gamma*self.__piLifetime)
        return t

    def Calculatet(self,s,s_final,t_i,v):
        t =  (np.abs(s) - np.abs(s_final))/v + t_i
        return t

    def Calculatez(self,s):
        return s

    def Calculatepz(self,p,px,py):
        pz = np.sqrt(p**2 - px**2 - py**2)
        return pz

    def CalculateE(self,p):
        E = np.sqrt(p**2 + nuSTORMTrfLine.__pimass**2)
        return E

    def Calculatev(self,p):
        E = self.CalculateE(p)
        beta = p/E
        v = beta * nuSTORMTrfLine.__sol
        return v

    def CalculateDcyPt(self,s,s_final,t_i,v,gamma):
        #t_f = time at entrance to production straight
        t_f = self.Calculatet(s=s,s_final=s_final,t_i=t_i,v=v)
        t_dcy = self.GenerateDcyTime(gamma)
        if t_dcy < (t_f - t_i):
            s_dcy = v*t_dcy - 50.
        else:
            s_dcy = -100.0
        return s_dcy

    def GeneratePion(self,**kwargs):
        eH = kwargs.get('eventHist')
        p0 = kwargs.get('p0')
        x  = kwargs.get('x')
        y  = kwargs.get('y')
        px = kwargs.get('px')
        py = kwargs.get('py')
        pz = kwargs.get('pz')
        t  = kwargs.get('t')
        weight = kwargs.get('weight')
        runNum = kwargs.get('runNum')
        eventNum = kwargs.get('eventNum')
        s = -50.
        s_final = 0.

        ## 1 - Generating pion at target

        if (eH == None or weight == None or runNum == None or eventNum == None or (p0 == None and (x == None or y == None or px == None or py == None or pz == None or t == None))):
            print("Input to function not valid!")
            return None
        elif p0 == None:
            Ppi = np.sqrt(px**2 + py**2 + pz**2)
            v = self.Calculatev(Ppi)
            z = self.Calculatez(s)
        else:
            Ppi = self.GeneratePiMmtm(p0)
            x,y,xp,yp = self.GenerateTrans(s)
            z = self.Calculatez(s)
            px = xp*Ppi
            py = yp*Ppi
            pz = self.Calculatepz(Ppi,px,py)
            t = self.GenerateTime()
            v = self.Calculatev(Ppi)

        pdg = self.__piPDG
        pi = particle.particle(runNum,eventNum,s,x,y,z,px,py,pz,t,weight,pdg)
        eH.addParticle("target",pi)


        ## 2 - Generating pion at production straight entrance and if pion decays in transfer line, also particles from pion decay at the decay point

        E = self.CalculateE(Ppi)
        gamma = E / self.__pimass
        s_dcy = self.CalculateDcyPt(s=s,s_final=s_final,t_i=t,v=v,gamma=gamma)

        if s_dcy == -100.:
            z_f = self.Calculatez(s_final)
            t_f = self.Calculatet(s,s_final,t,v)
            pi_prd = particle.particle(runNum,eventNum,s_final,x,y,z_f,px,py,pz,t_f,weight,pdg)
        else:
            z_dcy = self.Calculatez(s_dcy)
            t_dcy = self.Calculatet(s,s_dcy,t,v)
            DcyCoord = np.array([s_dcy,x,y,z_dcy])
            Dcy = PionDecay.PionDecay(ppi=Ppi)
            PEI = PionEventInstance.PionEventInstance(Ppi)
            P_mu, P_numu = PEI.Boost2nuSTORM(Dcy,Ppi)
            pi_dcy = particle.particle(runNum,eventNum,s_dcy,x,y,z_dcy,px,py,pz,t_dcy,weight,pdg)
            muon = self.GenerateMuon(runNum,eventNum,DcyCoord,P_mu,t_dcy,weight)
            numu = self.GenerateNu(runNum,eventNum,DcyCoord,P_numu,t_dcy,weight)
            eH.addParticle("pionDecay",pi_dcy)
            eH.addParticle("muonProduction",muon)
            eH.addParticle("piFlashNu",numu)
            pi_prd = particle.particle(-1,0,0.,0.,0.,0.,0.,0.,1.,0.,0.,pdg)

        eH.addParticle("productionStraight",pi_prd)
        return eH

    def GenerateMuon(self,runNum,eventNum,DcyCoord,P_mu,t,weight):
        s = DcyCoord[0]
        x = DcyCoord[1]
        y = DcyCoord[2]
        z = DcyCoord[3]
        px = P_mu[1][0]
        py = P_mu[1][1]
        pz = P_mu[1][2]
        pdg = self.__muPDG
        muon = particle.particle(runNum,eventNum,s,x,y,z,px,py,pz,t,weight,pdg)
        return muon

    def GenerateNu(self,runNum,eventNum,DcyCoord,P_numu,t,weight):
        s = DcyCoord[0]
        x = DcyCoord[1]
        y = DcyCoord[2]
        z = DcyCoord[3]
        px = P_numu[1][0]
        py = P_numu[1][1]
        pz = P_numu[1][2]
        pdg = -14
        numu = particle.particle(runNum,eventNum,s,x,y,z,px,py,pz,t,weight,pdg)
        return numu

    def Boost2nuSTORM(self, PPi, P_mu, P_numu):
        ''' Present approximation is muon propagates along z axis, so, boost only
            in preparation for later, include rotation matrix to transform from
            nustorm frame to frame with z axis along muon momentum and back '''

        EPi = np.sqrt(PPi**2 + self.__pimass**2)
        beta   = PPi / EPi
        gamma  = EPi / self.__pimass

        R    = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        Rinv = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])

        P_mu    = self.RotnBoost(P_mu, R, Rinv, gamma, beta)
        P_numu    = self.RotnBoost(P_numu, R, Rinv, gamma, beta)

        return P_mu, P_numu

    def RotnBoost(self, P, R, Rinv, gamma, beta):

        #print ("P[1] is ", P[1])
        #print ("P is ", P)
        p3 = np.dot(R, P[1])
        #print("    p3 is ", p3)

        Ec = P[0]
        Pc = p3[2]
        #print("    Ec is ", Ec, "   Pc is ", Pc)

        Ef = gamma * (Ec + beta * Pc)
        Pf = gamma * (Pc + beta * Ec)
        #print("    Ef is ", Ef, "   Pf is ", Pf)

        p3[2] = Pf
        p3    = np.dot(Rinv, p3)
        #print("    p3 is ", p3)

        Po    = [0., np.array([0., 0., 0.])]
        Po[0] = Ef
        Po[1] = p3
        #print("    Po is ", Po)

        return Po

    def BoostBack2RestFrame(self, PPi, P_mu, P_numu):
        ''' Present approximation is muon propagates along z axis, so, boost only
            in preparation for later, include rotation matrix to transform from
            frame with z axis along muon momentum to nustorm frame and back '''

        EPi = np.sqrt(PPi**2 + self.__pimass**2)
        beta   = PPi / EPi
        #print("beta: ",beta)
        gamma  = EPi / self.__pimass
        #print("gamma: ",gamma)

        R    = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        Rinv = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])

        P_mu    = self.NegRotnBoost(P_mu, R, Rinv, gamma, beta)
        P_numu    = self.NegRotnBoost(P_numu, R, Rinv, gamma, beta)

        return P_mu, P_numu

    def NegRotnBoost(self, P, R, Rinv, gamma, beta):

        #print ("P[1] is ", P[1])
        #print ("P is ", P)
        p3 = np.dot(R, P[1])
        #print("    p3 is ", p3)

        Ec = P[0]
        Pc = p3[2]
        #print("    Ec is ", Ec, "   Pc is ", Pc)

        Ef = gamma * (Ec - beta * Pc)
        Pf = gamma * (Pc - beta * Ec)
        #print("    Ef is ", Ef, "   Pf is ", Pf)

        p3[2] = Pf
        p3    = np.dot(Rinv, p3)
        #print("    p3 is ", p3)

        Po    = [0., np.array([0., 0., 0.])]
        Po[0] = Ef
        Po[1] = p3
        #print("    Po is ", Po)

        return Po

#--------  I/o methods:
    def GetTrfLineParams(_filename):
        TrfLineParams = pnds.read_csv(_filename)
        return TrfLineParams

#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)

    def CdVrsn(self):
        return 1.0

    def TrfLineLen(self):
        return deepcopy(self._TrfLineLen)

    def piAcc(self):
        return deepcopy(self._piAcc)

    def epsilon(self):
        return deepcopy(self._epsilon)

    def beta(self):
        return deepcopy(self._beta)

    def delT(self):
        return deepcopy(self._delT)

    def printParams(self):
        print(self._TrfLineParams)
