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
  _pAcc            = Momentum acceptance (%)
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
      pAcc         : Get momentum acceptance (%)
      epsilon      : Get acceptance/emittance (pi mm rad)
      beta         : Get beta function (mm)
      delT         : Get bunch length (ns)


  Generation, calcution methods and utilities:
      GenerateMmtm : Generates momentum centred on p0 (input float).
                     Parabolic distribution generated.
      GenerateTime : Generates time to enter transfer line; first particle leaves at t=0.
                     Uniform distribution generated.
      Calculatez   : Calculate z of decay given s, p and transfer line
                     and production-straight parameters.
                     **This version: z=s** (m)
      Calculatepz  : Calculate pz from px and py.
      GenerateTrans: Generate transverse phase space (x, y, xp, yp) given
                     representative emittance and beta.  Parabolic distributions
                     generated.
      GeneratePion : Generate a pion instance
      getTrfLineParams: Read csv file and generate Pandas dataframe
      printParams  : Print all paramters


Created on Tu 31Aug21. Version history:
----------------------------------------
 1.0: 31Aug21: First implementation

@author: MarvinPfaff
based on nuSTORMPrdStrght class
"""

import pandas as pnds
import numpy as np
from copy import deepcopy
import PionConst as PionConst
import MuonConst as MuonConst
import PionDecay as PionDecay
import particle as particle
import pion as pion
import muon as muon
import neutrino as neutrino
import Simulation as Simu

piCnst = PionConst.PionConst()
muCnst = MuonConst.MuonConst()


class nuSTORMTrfLine(object):
    __instance = None
    __pimass = piCnst.mass()/1000.
    __mumass = muCnst.mass()/1000.
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
            cls._pAcc            = cls._TrfLineParams.iat[1,2] / 100.
            cls._epsilon         = cls._TrfLineParams.iat[2,2]
            cls._beta            = cls._TrfLineParams.iat[3,2]
            cls._delT            = cls._TrfLineParams.iat[4,2]
            print('creation of new object: okay')
        return cls.__instance

    def __repr__(self):
        return "nuSTORMTrfLine()"

    def __str__(self):
        return "nuSTORMTrfLine: version=%g, length of transfer line=%g m," \
               "momentum acceptance=%g%%, transverse acceptance=%g pi mm rad, beta=%g mm, bunch length=%g ns. \n" % \
               (self.CdVrsn(), self.TrfLineLen(), self.pAcc(), self.epsilon(), self.beta(), self.delT() )

#--------  Simulation methods:
    def GenerateMmtm(self,p0):
        #import Simulation as Simu
        p = -99.
        dp = p0 * self._pAcc
        p  = p0 + Simu.getParabolic(dp)
        return p

    def GenerateTime(self):
        t = Simu.getRandom() * self._delT
        return t

    def GenerateTrans(self,s):
        #import Simulation as Simu
        r  = np.sqrt(self._epsilon*self._beta) / 1000.
        rp = np.sqrt(self._epsilon/self._beta)
        x  = Simu.getParabolic(r)
        y  = Simu.getParabolic(r)
        xp = Simu.getParabolic(rp)
        yp = Simu.getParabolic(rp)
        return x, y, xp, yp

    def Calculatet(self,s,s_final,v,t):
        t =  (np.abs(s) - np.abs(s_final))/v + t
        return t

    def Calculatez(self,s):
        return s

    def Calculatepz(self,p,px,py):
        pz = np.sqrt(p**2 - px**2 - py**2)
        return pz

    def Calculatev(self,p):
        E = np.sqrt(p**2 + nuSTORMTrfLine.__pimass**2)
        beta = p/E
        v = beta * nuSTORMTrfLine.__sol
        return v

    def CalculateDcyPt(self,s,s_final,t,v):
        t1 = self.Calculatet(s=s,s_final=s_final,t=t,v=v)    #time at entrance to production straight
        Dcy = PionDecay.PionDecay()
        t_decay = Dcy.getLifetime()
        if t_decay < (t1 - t):
            s_decay = v*t_decay
        else:
            s_decay = -100.0
        return t_decay

    def GeneratePion(self,**kwargs):
        p0 = kwargs.get('p0')
        s  = kwargs.get('s')
        x  = kwargs.get('x')
        y  = kwargs.get('y')
        z  = kwargs.get('z')
        px = kwargs.get('px')
        py = kwargs.get('py')
        pz = kwargs.get('pz')
        t  = kwargs.get('t')
        weight = kwargs.get('weight')
        runNum = kwargs.get('runNum')
        eventNum = kwargs.get('eventNum')
        s_final = kwargs.get('s_final')
        s_decay = kwargs.get('s_decay')
        Dcy = kwargs.get('Dcy')

        #Adjust if statement --> also add weight to be mandatory!!!
        #if (s_final == None || weight == None || (p0 == None && (runNum == None || eventNum == None || s == None || x == None || y == None || z == None || px == None || py == None || pz == None || t == None || weight == None) && s_decay == None)):
        #    print("Input to function not valid!")
        #    return -200.0, pion.pion(0,0,0.,0.,0.,0.,0.,0.,0.,0.,0.)

        #elif p0 == None:
        if p0 == None:
            p = np.sqrt(px**2 + py**2 + pz**2)
            v = self.Calculatev(p)
            if s_final != s:
                t = self.Calculatet(s=s,s_final=s_final,t=t,v=v)
                z = self.Calculatez(s)

                if s_decay == s_final:
                    runNum = -1
                    #Input rest here!!

        else:
            s = s_final
            s_final = 0.
            p = self.GenerateMmtm(p0)
            x,y,xp,yp = self.GenerateTrans(s)
            z = self.Calculatez(s)
            px = xp*p
            py = yp*p
            pz = self.Calculatepz(p,px,py)
            t = self.GenerateTime()
            v = self.Calculatev(p)

            s_decay = self.CalculateDcyPt(s=s,s_final=s_final,t=t,v=v)

        pdg = 211
        mass = nuSTORMTrfLine.__pimass
        pi = particle.particle(runNum,eventNum,s,x,y,z,px,py,pz,t,weight,mass,pdg)
        return s_decay, pi

    def GenerateMuon(self,runNum,eventNum,DcyCoord,P_mu,t,weight):
        s = DcyCoord[0]
        x = DcyCoord[1]
        y = DcyCoord[2]
        z = DcyCoord[3]
        px = P_mu[1]
        py = P_mu[2]
        pz = P_mu[3]
        mass = MuonConst.mass()
        pdg =
        muon = particle.particle(runNum,eventNum,s,x,y,z,px,py,pz,t,weight,mass,pdg)
        return muon

    def GenerateNu(self,runNum,eventNum,DcyCoord,P_numu,t,weight):
        s = DcyCoord[0]
        x = DcyCoord[1]
        y = DcyCoord[2]
        z = DcyCoord[3]
        px = P_numu[1]
        py = P_numu[2]
        pz = P_numu[3]
        mass =
        pdg =
        numu = particle.particle(runNum,eventNum,s,x,y,z,px,py,pz,t,weight,mass,pdg)
        return numu

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

    def pAcc(self):
        return deepcopy(self._pAcc)

    def epsilon(self):
        return deepcopy(self._epsilon)

    def beta(self):
        return deepcopy(self._beta)

    def delT(self):
        return deepcopy(self._delT)

    def printParams(self):
        print(self._TrfLineParams)
