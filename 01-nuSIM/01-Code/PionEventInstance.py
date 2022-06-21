#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class PionEventInstance:
============================

  Generates single muon production event along the nuSTORM production
  straight.

  Class attributes:
  -----------------
  __PionDecay: pion decay class
  --np       : numpy class

  Instance attributes:
  --------------------
  _ppi          : Pion momentum: i/p argument at instance creation
  _pmu          : Muon momentum;
  _TrcSpcCrd    : Trace space (s, x, y, z, x', y') in numpy array at
                  point of decay
  _LclTrcSpcCrd : (Local) Trace space (s, x, y, z, x', y') in numpy
                  array at point where pion is generated
  _ppiGen       : Generated pion momentum
  _P_mu         : Muon 4 momentum: (E, array(px, py, pz)), GeV
  _P_numu       : Muon-neutrino 4 momentum: (E, array(px, py, pz)), GeV

  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates decay instance
      __repr__ : One liner with call.
      __str__  : Dump of values of decay

  Get/set methods:
    getpppi              : Returns pion momentum as real
    getpmu               : Returns muon momentum as real
    getTraceSpaceCoord   : Returns trace space: (s, x, y, z, x', y') (m) at pion decay
    getLclTraceSpaceCoord: Returns (local) trace space: (s, x, y, z, x', y') (m) at pion generation
    getppiGen            : Returns generated pion momentum (GeV)
    getmu4mmtm           : Returns muon 4 momentum: (E, array(px, py, pz)), GeV
    getnumu4mmtm         : Returns muon-neutrino 4 momentum: (E, array(px, py, pz)), GeV

  General methods:
    CreateMuon           : Manager for pion decay, returns z (m) of decay (real), P_mu, P_numu, RestFrame
                           Restframe contains a dump of the instance attributes of the PionDecay class
    GenerateDcyPhaseSpace: Generates trace-space position of decay
    GenerateLongiPos     : Returns s, z of decay
    Boost2nuSTORM        : Boots to nuSTORM rest frame -- i.e. boost to ppi
    RotnBoost            : Operator; rotates and boosts rest-frame coordinates to nuSTORM frame

Created on Tue 30Mar21;02:26: Version history:
----------------------------------------------
 1.0: 30Mar21: First implementation - checked on 12 May that the decay length for a 6 GeV pion beam
 is as expected, by fitting the distribution of decay lengths

 1.1: 20Oct21: Correct problem created by adding different momentum distributions for pions and muons
@author: kennethlong
@author: PaulKyberd

 1.2: 23Nov21: Pion decay cut on production straightLength causes normalisation problems. Modify to take all decays
@author: PaulKyberd

 1.3: 20Jan22: Include possibility to take original pion phase space from pion at target in event history and exploit
               TLorentzVector class from PyROOT to do bost to nuSTORM frame
 @author: MarvinPfaff

 1.3: 06Jun22: Output two trace spaces, one at pion generation (Coord, _LclTrcSpcCrd) and one at pion decay
               (DcyCoord, _TrcSpcCrd) and get nuSTORM constants from nuSTORMConst.py
 @author: MarvinPfaff
"""

from copy import deepcopy
from ROOT import TLorentzVector as T4V
import nuSTORMPrdStrght as nuPrdStrt
#import nuSTORMTrfLineCmplx as nuTrf
import nuSTORMConst
import PionDecay as PionDecay
import MuonConst as MuonConst
import PionConst as PionConst
import numpy as np
import math
import os

muCnst = MuonConst.MuonConst()
piCnst = PionConst.PionConst()
nuSTRMCnst = nuSTORMConst.nuSTORMConst()
nuStrt = nuPrdStrt.nuSTORMPrdStrght('11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
nuSIMPATH = os.getenv('nuSIMPATH')
#trfCmplxFile = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-TrfLineCmplx-Params-v1.0.csv')
#nuTrLnCmplx = nuTrf.nuSTORMTrfLineCmplx(trfCmplxFile)
tlCmplxLength = nuSTRMCnst.TrfLineCmplxLen()
tlCmplxAngle  = nuSTRMCnst.TrfLineCmplxAng()

class PionEventInstance:

    __mumass = muCnst.mass()/1000.
    __pimass = piCnst.mass()/1000.
    __sol    = muCnst.SoL()

    __Debug  = False

#--------  "Built-in methods":
    def __init__(self, ppi=6., **kwargs):

        self._pion = kwargs.get('particleTar')
        self._pionLcl = kwargs.get('particleTarLocal')
        if self._pion != None:
            #self._mode marks whether phase space should be taken from input or created according to a distribution
            self._mode = 'input'
        elif self._pionLcl != None:
            self._mode = 'input local'
            self._pion = self._pionLcl
        else:
            self._ppi = ppi
            self._mode = 'random'
        self._phi=-10.0
        self._costheta=-10.0
        self._lifetime = 0.0
        self._P_pi=np.array([0., np.array([0., 0., 0.])])
        self._tlAngle = tlCmplxAngle * math.pi / 180.
        self._TrcSpcCrd, self._LclTrcSpcCrd, self._ppiGen, self._P_mu, self._P_numu = self.CreateMuon()

        return

    def __repr__(self):
        return "PionEventInstance(ppi)"

    def __str__(self):
        return "PionEventInstance: p_pi (GeV) = %g, s (m) = %g, z (m) = %g, generated momentum=%g, \r\n \
                P_pi (%g, [%g, %g, %g]), \r\n \
                P_mu (%g, [%g, %g, %g]), \r\n \
                P_numu (%g, [%g, %g, %g]), \r\n" % \
            (self._ppi, self._TrcSpcCrd[0], self._TrcSpcCrd[3], self._ppiGen, \
             self._P_pi[0], self._P_pi[1][0],self._P_pi[1][1],self._P_pi[1][2], \
             self._P_mu[0], self._P_mu[1][0],self._P_mu[1][1],self._P_mu[1][2], \
             self._P_numu[0], self._P_numu[1][0],self._P_numu[1][1],self._P_numu[1][2] )

#--------  Generation of neutrino-creation event:
#.. Manager:
    def CreateMuon(self):
#        print ("Starting create muon")

        PrdStrghtLngth = nuSTRMCnst.ProdStrghtLen()
#        #.. Prepare--get muon decay instance in pion rest frame:
#        z = 2.* PrdStrghtLngth
#        Dcy = 0
#        if PionEventInstance.__Debug:
#            print("PionEventInstance.CreatePion: find valid decay")
#        while z > PrdStrghtLngth:
#            if isinstance(Dcy, PionDecay.PionDecay):
#                del Dcy
        if self._mode == 'random':
            Ppi0 = self.getppi()
#   Comment out this line for distribution
#            Ppi = Ppi0
#    Comment out this line and we get single energy
            Ppi = nuStrt.GeneratePiMmtm(Ppi0)
            Epi   = np.sqrt(Ppi**2 + PionEventInstance.__pimass**2)
        else:
            px = self._pion.p()[1][0]
            py = self._pion.p()[1][1]
            pz = self._pion.p()[1][2]
            Ppi = np.sqrt(px**2 + py**2 + pz**2)
            self._ppi = Ppi
            Epi = self._pion.p()[0]

        beta  = Ppi / Epi
        gamma = Epi / PionEventInstance.__pimass
        v    = beta * PionEventInstance.__sol
        if (self.__Debug): 
            print(f"PionEventInstance.CreatePion: Ppi {Ppi}   Epi {Epi}   beta {beta}   gamma {gamma}   v {v}")
            gammaChck = 1.0/math.sqrt(1-beta*beta)
            print(f"     delta gamma {gammaChck - gamma}\n\n")

# no TMax at present
#        Tmax = nuStrt.ProdStrghtLen() / (gamma * v)
#        Dcy = PionDecay.PionDecay(Tmax=Tmax)
        Dcy = PionDecay.PionDecay()
        self._phi = Dcy.getphi()
        self._costheta = Dcy.getcostheta()
        if (self.__Debug): print("PionEventInstance.CreatePion: Dcy ", Dcy)
        DcyCoord, Coord, DirCos = self.GenerateDcyPhaseSpace(Dcy,Ppi)
        if self._mode == 'random':
            self._P_pi[0] = Epi 
            self._P_pi[1][0] = Ppi*DirCos[0]
            self._P_pi[1][1] = Ppi*DirCos[1]
            self._P_pi[1][2] = Ppi*DirCos[2]
        self._lifetime = Dcy.getLifetime()
        if (self.__Debug): print("PionEventInstance.CreatePion: Ppi ", Ppi)
        if (self.__Debug): print(f"PionEventInstance.CreatePion: DcyCoord  {DcyCoord} \n          Coord {Coord} \n          DirCos {DirCos}")
        z = DcyCoord[3]
#  This test is no longer valid - may re-introduce
#        if z > PrdStrghtLngth:
#            print("PionEventInstance.CreatePion Alarm:", z)
        if PionEventInstance.__Debug:
            print("PionEventInstance.CreatePion: decay at z =", z)
            print("----> Dcy:", Dcy.__str__())

        #.. Boost to nuSTORM frame:
        if PionEventInstance.__Debug:
            print(f"PionEventInstance.CreateMuon: rotate and boost to nuSTORM rest frame:Ppi is {Ppi}")
        P_mu, P_numu = self.Boost2nuSTORM(Dcy,Ppi, DirCos)

        del Dcy

#        print ("Ending create muon")
        return DcyCoord, Coord, Ppi, P_mu, P_numu

#.. Trace space coordinate generation: array(s, x, y, z, x', y')
    def GenerateDcyPhaseSpace(self, Dcy, Ppi):
 #       print ("generate phase space start")
        dcyCoord = np.array([0., 0., 0., 0., 0., 0.])
        coord = np.array([0., 0., 0., 0., 0., 0.])
        #.. longitudinal position, "s", z:
#        dcyCoord[0] = self.GenerateLongiPos(Dcy, Ppi)
# using  s = (beta*c) * (gamma*t) = p/E * c * E/m * lifetime = p*c*lifetime/m
        dcyCoord[0] =Ppi*PionEventInstance.__sol*Dcy.getLifetime()/PionEventInstance.__pimass

        R, Rinv, BeamPos, theta = self.BeamDir(dcyCoord[0], Ppi)
# here we set the Beam Twiss parameters and Energy spread

        if self._mode == 'random':
            xl, yl, xpl, ypl = nuStrt.GenerateTrans(dcyCoord[0])
            s  = 0.
            zl = s - tlCmplxLength
        elif (self._mode == 'input' and self._pion.traceSpace().s() < tlCmplxLength):
            #transform the global coordinates of pion at target back to tlLocal to get the position and momentum spread with respect to BeamPos & BeamDir
            s = self._pion.traceSpace().s()
            xg = self._pion.x()
            yg = self._pion.y()
            zg = self._pion.z()
            pxg = self._pion.p()[1][0]
            pyg = self._pion.p()[1][1]
            pzg = self._pion.p()[1][2]
            xl, yl, zl, pxl, pyl, pzl = self.glbltoTl(xg, yg, zg, pxg, pyg, pzg)
            xpl = pxl/pzl
            ypl = pyl/pzl
        elif (self._mode == 'input local' or (self._mode == 'input' and self._pion.traceSpace().s() >= tlCmplxLength)):
            s = self._pion.traceSpace().s()
            xl = self._pion.x()
            yl = self._pion.y()
            zl = self._pion.z()
            pxl = self._pion.p()[1][0]
            pyl = self._pion.p()[1][1]
            pzl = self._pion.p()[1][2]
            xpl = pxl/pzl
            ypl = pyl/pzl

        coord[0] = s
        coord[1] = xl
        coord[2] = yl
        coord[3] = zl
        coord[4] = xpl
        coord[5] = ypl

        dcyCoord[1] = xl + BeamPos[0]
        dcyCoord[2] = yl + BeamPos[1]
        # [3] is the z position which is calculated from the decay distance and information about the
        # beam trajectory.
        dcyCoord[3] = BeamPos[2]
        dcyCoord[4] = xpl
        dcyCoord[5] = ypl

        p0    = np.array([0., 0., 0.])
        p1    = np.array([0., 0., 0.])
        p0[0] = xpl
        p0[1] = ypl
        p0[2] = math.sqrt(1. - xpl**2 - ypl**2)

        p1    = R.dot(p0)

        if PionEventInstance.__Debug:
            print("         ----> Direction cosines from trace space: ", p0)
            print("         ----> R: ",  R[0], R[1], R[2])
            print("         ----> Direction cosines rotated: ", p1)

        return dcyCoord, coord, p1
        #return dcyCoord, Ppi

#.. Beam position, direction and corresponding rotation operator:

    def BeamDir(self, s, Ppi):
        PrdStrghtLngth = nuSTRMCnst.ProdStrghtLen()
        Circumference  = nuSTRMCnst.Circumference()
        ArcLen         = nuSTRMCnst.ArcLen()

        #s_ring takes into account that s = 0. at target and s_ring=z for production straight
        s_ring         = s-tlCmplxLength
        #where corresponds to location within the ring
        if (s_ring <= 0.):
            where      = s_ring
        else:
            where      = s_ring%Circumference
        ArcRad         = ArcLen/math.pi

        if ( 0. >= where):
            if PionEventInstance.__Debug:
                print("         ---->", \
                      "BeamDir: the pion decays in the transfer line complex")

            theta         =  self._tlAngle          #Angle with respect to Z axis
            z             =  math.cos(theta)*where
            x             = -math.sin(theta)*where
            y             = 0.
            BeamPos       = [x, y, z]

        elif ( PrdStrghtLngth >= where > 0.):
            if PionEventInstance.__Debug:
                print("         ---->", \
                      "BeamDir: the pion decays in the production straight")

            theta         = 0.     #Angle with respect to Z axis
            BeamPos       = [0., 0., where]

        elif (PrdStrghtLngth+ArcLen >= where > PrdStrghtLngth):
            if PionEventInstance.__Debug:
                print("         ---->", \
                      "BeamDir: the pion decays in the first arc")

            ArcLenCovered = where - PrdStrghtLngth
            theta         = math.pi*ArcLenCovered/ArcLen
            BeamPos       = [ \
                              (ArcRad-ArcRad* math.cos(theta)),       \
                              0.,                                     \
                              PrdStrghtLngth+ArcRad*math.sin(theta)]

        elif (2*PrdStrghtLngth+ArcLen >= where > PrdStrghtLngth+ArcLen):
            if PionEventInstance.__Debug:
                print("         ---->", \
                      "BeamDir: the pion decays in the return straight")

            theta         = math.pi
            BeamPos       = [2.*ArcRad, \
                             0., \
                             PrdStrghtLngth-(where-ArcLen-PrdStrghtLngth)]

        elif (2*PrdStrghtLngth+2*ArcLen>=where>=2*PrdStrghtLngth+ArcLen):
            if PionEventInstance.__Debug:
                print("         ---->", \
                      "BeamDir: the pion decays in the second/return arc")

            ArcLenCovered = where - 2*PrdStrghtLngth - ArcLen
            theta         = math.pi + math.pi*ArcLenCovered/ArcLen
            BeamPos       = [ \
                              (ArcRad-ArcRad*math.cos(theta)), \
                              0., \
                              ArcRad*math.sin(theta)]

        # R=rotation with respect to y axis through theta angle
        R    = np.array([ \
                          [math.cos(theta),  0., math.sin(theta)], \
                          [0.,               1., 0.],              \
                          [-math.sin(theta), 0., math.cos(theta)]])
        Rinv = np.array([ \
                          [math.cos(theta),  0.,-math.sin(theta)], \
                          [0.,               1., 0.],              \
                          [math.sin(theta),  0., math.cos(theta)]])

        if PionEventInstance.__Debug:
            print('               ----> theta, s:', theta, s)

        return R, Rinv, BeamPos, theta

#.. Boost from muon rest frame to nuSTORM frame:
    def Boost2nuSTORM(self, Dcy, PPi, DirCos):
        ''' Present approximation is muon propagates along z axis, so, boost only
            in preparation for later, include rotation matrix to transform from
            nustorm frame to frame with z axis along muon momentum and back '''

#        print ("PPi on Start boost to nuStorm", PPi)


        EPi = np.sqrt(PPi**2 + PionEventInstance.__pimass**2)
        pB  = PPi*DirCos
        pi4v = T4V.TLorentzVector(pB[0], pB[1], pB[2], EPi)
        b    = pi4v.BoostVector()

        #beta   = PPi / EPi
        #gamma  = EPi / PionEventInstance.__pimass

        #R    = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        #Rinv = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        if PionEventInstance.__Debug:
            print("PionEventInstance.Boost2nuSTORM: boost parameters:")
            print("----> PPi, EPi, boost:", PPi, EPi, "; ", b[0], b[1], b[2])    # OK with P_mu=8

        # Treat decay components:
        P_mu    = Dcy.get4vmu()
        P_mu[0] = P_mu[0]/1000.
        P_mu[1] = P_mu[1]/1000.

        P_mu4v  = T4V.TLorentzVector(P_mu[1][0], \
                                    P_mu[1][1], \
                                    P_mu[1][2], \
                                    P_mu[0] )
        if PionEventInstance.__Debug:
            print("            ----> Rest frame P_mu4v (GeV):",
                  P_mu4v.E(), P_mu4v.Px(), P_mu4v.Py(), P_mu4v.Pz())

        P_mu4v.Boost(b)

        if PionEventInstance.__Debug:
            print("            ----> nuSTORM frame P_mu4v (GeV):",
                  P_mu4v.E(), P_mu4v.Px(), P_mu4v.Py(), P_mu4v.Pz())

        P_numu    = Dcy.get4vnumu()
        P_numu[0] = P_numu[0]/1000.
        P_numu[1] = P_numu[1]/1000.
        P_numu4v  = T4V.TLorentzVector(P_numu[1][0], \
                                       P_numu[1][1], \
                                       P_numu[1][2], \
                                       P_numu[0] )

        if PionEventInstance.__Debug:
            print("            ----> Rest frame P_numu4v (GeV):",
                  P_numu4v.E(), P_numu4v.Px(), P_numu4v.Py(), P_numu4v.Pz())

        P_numu4v.Boost(b)

        if PionEventInstance.__Debug:
            print("            ----> nuSTORM frame P_numu4v (GeV):",
                  P_numu4v.E(), P_numu4v.Px(), P_numu4v.Py(), P_numu4v.Pz())

        P_mu[0]    = P_mu4v.E()
        P_mu[1][0] = P_mu4v.Px()
        P_mu[1][1] = P_mu4v.Py()
        P_mu[1][2] = P_mu4v.Pz()

        P_numu[0]    = P_numu4v.E()
        P_numu[1][0] = P_numu4v.Px()
        P_numu[1][1] = P_numu4v.Py()
        P_numu[1][2] = P_numu4v.Pz()

 #       print ("End boost to nuStorm")

        return P_mu, P_numu


# transform x,y,z and px,py,pz co-ordinates from the global co-ordinates to the
# transfer line local ones

    def glbltoTl(self, xg, yg, zg, pxg, pyg, pzg):


        xl = xg + zg*math.tan(self._tlAngle)
        yl = yg
        zl = zg/math.cos(self._tlAngle)

        pxl = pxg*math.cos(self._tlAngle) - pzg*math.sin(self._tlAngle)
        pyl = pyg
        pzl = pxg*math.sin(self._tlAngle) + pzg*math.cos(self._tlAngle)

        return xl, yl, zl, pxl, pyl, pzl


#--------  get/set methods:
    def getLifetime(self):
        return deepcopy(self._lifetime)

    def getppi(self):
        return deepcopy(self._ppi)

    def getphi(self):
        return deepcopy(self._phi)

    def getcostheta(self):
        return deepcopy(self._costheta)

    def getTraceSpaceCoord(self):
        return deepcopy(self._TrcSpcCrd)

    def getLclTraceSpaceCoord(self):
        return deepcopy(self._LclTrcSpcCrd)

    def getpi4mmtm(self):
        return deepcopy(self._P_pi)

    def getmu4mmtm(self):
        return deepcopy(self._P_mu)

    def getnumu4mmtm(self):
        return deepcopy(self._P_numu)

    def getppiGen(self):
        return deepcopy(self._ppiGen)
