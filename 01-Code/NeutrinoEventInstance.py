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

      
  Instance attributes:
  --------------------
  _pmu      : Muon momentum; i/p argument at instance creation
  _ct       : Decay time * speed-of-light
  _TrcSpcCrd: Trace space (s, x, y, z, x', y') in numpy array at 
              point of decay
  _pmuGen   : Generated muon momentum
  _pmuDirCos: np.array([d, d, d]) with three direction cosines
  _P_e      : Electron 4 momentum: (E, array(px, py, pz)), GeV
  _P_nue    : Electron-neutrino 4 momentum: (E, array(px, py, pz)), GeV
  _P_numu   : Muon-neuutrino 4 momentum: (E, array(px, py, pz)), GeV

      
  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates decay instance; calls CreateNeutrinos method to do it
      __repr__ : One liner with call
      __str__  : Dump of values of decay and beam params at decay point

  Get/set methods:
    getpmu            : Returns nominal muon momentum (GeV)
    getTraceSpaceCoord: Returns trace space: (s, x, y, z, x', y') (m)
    getpmuGen         : Returns generated muon momentum (GeV)
    getPb             : Returns three vector (np.array) of beam momentum (GeV)
    gete4mmtm         : Returns electron 4 momentum: (E, array(px, py, pz))
                            (GeV)
    getnue4mmtm       : Returns electron-neutrino 4 momentum: 
                            (E, array(px, py, pz)), GeV
    getnumu4mmtm      : Returns muon-neutrino 4 momentum: 
                            (E, array(px, py, pz)), GeV


  General methods:
    CreateNeutrinos      : Manager for neutrino decay.  Returns:
                           ct, [trace space], pmuGen, [DirCos], [P_e],
                           [P_nue], [P_numu]
    GenerateDcyPhaseSpace: Generates trace-space position of decay.  Returns:
                           [trace space], [DirCos]
    GenerateLongiPos     : Returns s of decay
    BeamDir              : Returns position and momentum of the beam,
                           rotation operator corresponding to the direction 
                           of beam momentum and angle of beam wrt z axis at
                           decay.  Returns:
                           [R[3][3]], [Rinv[3][3]], [Pos. of decay], theta
    Boost2nuSTORM        : Boost to nuSTORM rest frame -- i.e. boost to pmuGen.
                           Uses TLorentVector class from ROOT.  Returns:
                           [P_e], [P_nue], [P_numu]


Created on Sat 16Jan21;02:26: Version history:
----------------------------------------------
 1.3: 28Aug21: KL: Review and tidy Omar's code and exploit TLorentzVector
               class from PyROOT to do bost to nuSTORM frame
 1.2: 18Jun21: Pass nuSTORM production straight parameter file name
 1.1: 03Apr21: Fix error in relativistic treatment muon lifetime
 1.0: 16Jan21: First implementation

@author: kennethlong
"""

from copy import deepcopy

from ROOT import TLorentzVector as T4V

import nuSTORMPrdStrght as nuPrdStrt
import MuonDecay as MuonDecay
import MuonConst as MuonConst
import numpy as np
import math
import PionEventInstance as piEvt
import Simulation as Simu


muCnst = MuonConst.MuonConst()
Ppi0=8
piEvt  = piEvt.PionEventInstance(Ppi0)

    
class NeutrinoEventInstance:

    __mumass = muCnst.mass()/1000.
    __sol    = muCnst.SoL()

    __Debug  = False

#--------  "Built-in methods":
    def __init__(self, piEvt, filename=None):

        nuStrt = nuPrdStrt.nuSTORMPrdStrght(filename)
  

        self._pmu = pmu
        self._ct, self._TrcSpcCrd, self._pmuGen, self._pmuDirCos,  \
            self._P_e, self._P_nue, self._P_numu \
                = self.CreateNeutrinos(nuStrt, PiEvt)

        return

    def __repr__(self):
        return "NeutrinoEventInstance(pmu)"

    def __str__(self):
        return "NeutrinoEventInstance: p_mu (GeV) = %g, " \
               " s (m) = %g, z (m) = %g,\n" \
               " Generated momentum=%g, direction cosines=[%g, %g, %g]\n" \
               " P_e (%g, [%g,, %g, %g]), \n" \
               " P_nue (%g, [%g,, %g, %g]), \n" \
               " P_numu (%g, [%g,, %g, %g])" % \
       (self._pmu, self._TrcSpcCrd[0], self._TrcSpcCrd[3], self._pmuGen, \
        self._pmuDirCos[0], self._pmuDirCos[1], self._pmuDirCos[2], \
       self._P_e[0], self._P_e[1][0], self._P_e[1][1], self._P_e[1][2], \
       self._P_nue[0], self._P_nue[1][0],self._P_nue[1][1],self._P_nue[1][2], \
       self._P_numu[0], self._P_numu[1][0],self._P_numu[1][1],self._P_numu[1][2] )
    
#--------  Generation of neutrino-creation event:
#.. Manager:
    def CreateNeutrinos(self, nuStrt, piEvt):
        PrdStrghtLngth = nuStrt.ProdStrghtLen()
        Circumference  = nuStrt.Circumference()
        ArcLen         = nuStrt.ArcLen()       
        ArcRad         = ArcLen / math.pi
        
        tpi = piEvt.gettpi()
        spi = PiEvt.getTraceSpaceCoord()[0]
        Emu = piEvt.getmu4mmtm()[0]
        Pmu = piEvt.getmu4mmtm()[1]
        

        if NeutrinoEventInstance.__Debug:
           print(" NeutrinoEventInstance.CreateNeutrinos: entered")

        #.. Prepare--get neutrino decay instance in muon rest frame:
        z = 2.* PrdStrghtLngth
        Dcy = 0
        if NeutrinoEventInstance.__Debug:
            print("    ----> Find valid decay")
            
        if isinstance(Dcy, MuonDecay.MuonDecay):
            del Dcy
        Dcy = MuonDecay.MuonDecay()

        Pmu0             = self.getpmu()
        PmuGen           = nuStrt.GenerateMmtm(Pmu0)
        DcyCoord, DirCos = self.GenerateDcyPhaseSpace(Dcy, PmuGen, nuStrt)

        z  = DcyCoord[3]
        ct = Dcy.getLifetime()
        s  = DcyCoord[0]

        if z > (PrdStrghtLngth+ArcRad+1.):
            print("     ----> !!!! CreateNeutrinos Alarm:", z)

        if NeutrinoEventInstance.__Debug:
            print("    ----> Decay at z =", z)
            
        if NeutrinoEventInstance.__Debug:
            print("    ----> Rotate and boost to nuSTORM rest frame")

        P_e, P_nue, P_numu = self.Boost2nuSTORM(Dcy, PmuGen, DirCos)

        if NeutrinoEventInstance.__Debug:
            print("     ----> Decay products in nuSTORM frame:")
            print("         ----> P_e   :", P_e)
            print("         ----> P_nue :", P_nue)
            print("         ----> P_numu:", P_numu)

        del Dcy
        
        return ct, DcyCoord, PmuGen, DirCos, P_e, P_nue, P_numu

#.. Trace space coordinate generation: array(s, x, y, z, x', y')
    def GenerateDcyPhaseSpace(self, Dcy, Pmu, nuStrt):

        coord = np.array([0., 0., 0., 0., 0., 0.])

        #.. longitudinal position, "s", z:
        coord[0] = self.GenerateLongiPos(Dcy, Pmu)

        R, Rinv, BeamPos, theta = self.BeamDir(coord[0], nuStrt, Pmu)
        
        x, y, xp, yp = nuStrt.GenerateTrans(coord[0])
        
        coord[1] = x + BeamPos[0]
       
        coord[2] = y + BeamPos[1]
  
        coord[3] = BeamPos[2]
      
        coord[4] = xp
        coord[5] = yp

        p0    = np.array([0., 0., 0.])
        p1    = np.array([0., 0., 0.])
        p0[0] = xp
        p0[1] = yp
        p0[2] = math.sqrt(1. - xp**2 - yp**2)

        p1    = R.dot(p0)

        if NeutrinoEventInstance.__Debug:
            print("         ----> Direction cosines from trace space: ", p0)
            print("         ----> R: ",  R[0], R[1], R[2])
            print("         ----> Direction cosines rotated: ", p1)
        
        return coord, p1

#.. Trace space coordinate generation:
    def GenerateLongiPos(self, Dcy, Pmu):
        Emu   = np.sqrt(Pmu**2 + NeutrinoEventInstance.__mumass**2)
        beta  = Pmu / Emu
        gamma = Emu / NeutrinoEventInstance.__mumass
        
        v    = beta * NeutrinoEventInstance.__sol

        s   = v * gamma * Dcy.getLifetime()
        
        return s

#.. Beam position, direction and corresponding rotation operator:

    def BeamDir(self, s, nuStrt, Pmu):
        PrdStrghtLngth = nuStrt.ProdStrghtLen()
        Circumference  = nuStrt.Circumference()
        ArcLen         = nuStrt.ArcLen()

        where          = s%Circumference
        ArcRad         = ArcLen/math.pi
        
        if ( PrdStrghtLngth>=where):
            if NeutrinoEventInstance.__Debug:
                print("         ---->", \
                      "BeamDir: the muon decays in the production straight")

            theta         = 0.     #Angle with respect to Z axis
            BeamPos       = [0., 0., where]

        elif (PrdStrghtLngth+ArcLen >= where > PrdStrghtLngth):
            if NeutrinoEventInstance.__Debug:
                print("         ---->", \
                      "BeamDir: the muon decays in the first arc")

            ArcLenCovered = where - PrdStrghtLngth
            theta         = math.pi*ArcLenCovered/ArcLen
            BeamPos       = [ \
                              (ArcRad-ArcRad* math.cos(theta)),       \
                              0.,                                     \
                              PrdStrghtLngth+ArcRad*math.sin(theta)]

        elif (2*PrdStrghtLngth+ArcLen >= where > PrdStrghtLngth+ArcLen):
            if NeutrinoEventInstance.__Debug:
                print("         ---->", \
                      "BeamDir: the muon decays in the return straight")

            theta         = math.pi
            BeamPos       = [2.*ArcRad, \
                             0., \
                             PrdStrghtLngth-(where-ArcLen-PrdStrghtLngth)]

        elif (2*PrdStrghtLngth+2*ArcLen>=where>=2*PrdStrghtLngth+ArcLen):
            if NeutrinoEventInstance.__Debug:
                print("         ---->", \
                      "BeamDir: the muon decays in the second/return arc")

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
                        
        if NeutrinoEventInstance.__Debug: 
            print('               ----> theta, s:', theta, s)
                
        return R, Rinv, BeamPos, theta


#.. Boost from muon rest frame to nuSTORM frame:
    def Boost2nuSTORM(self, Dcy, Pmu, DirCos):
       
        Emu = np.sqrt(Pmu**2 + NeutrinoEventInstance.__mumass**2)
        pB  = Pmu*DirCos
        mu4v = T4V.TLorentzVector(pB[0], pB[1], pB[2], Emu)
        b    = mu4v.BoostVector()
        
        
        if NeutrinoEventInstance.__Debug:
            print("        ----> Boost2nuSTORM: boost parameters:")
            print("            ----> Pmu, Emu, boost:", Pmu, Emu, \
                  "; ", b[0], b[1], b[2])

        # Treat decay components:
        P_e    = Dcy.get4ve()
        P_e[0] = P_e[0]/1000.
        P_e[1] = P_e[1]/1000.
     
        P_e4v  = T4V.TLorentzVector(P_e[1][0], \
                                    P_e[1][1], \
                                    P_e[1][2], \
                                    P_e[0] )
        if NeutrinoEventInstance.__Debug:
            print("            ----> Rest frame P_e4v (GeV):",
                  P_e4v.E(), P_e4v.Px(), P_e4v.Py(), P_e4v.Pz())

        P_e4v.Boost(b)

        if NeutrinoEventInstance.__Debug:
            print("            ----> nuSTORM frame P_e4v (GeV):", 
                  P_e4v.E(), P_e4v.Px(), P_e4v.Py(), P_e4v.Pz())
        
        P_nue     = Dcy.get4vnue()
        P_nue[0]  = P_nue[0]/1000.
        P_nue[1]  = P_nue[1]/1000.
    
        P_nue4v   = T4V.TLorentzVector(P_nue[1][0], \
                                       P_nue[1][1], \
                                       P_nue[1][2], \
                                       P_nue[0] )

        if NeutrinoEventInstance.__Debug:
            print("            ----> Rest frame P_nue4v (GeV):",
                  P_nue4v.E(), P_nue4v.Px(), P_nue4v.Py(), P_nue4v.Pz())

        P_nue4v.Boost(b)
        
        if NeutrinoEventInstance.__Debug:
            print("            ----> nuSTORM frame P_nue4v (GeV):", 
                  P_nue4v.E(), P_nue4v.Px(), P_nue4v.Py(), P_nue4v.Pz())
        
        P_numu    = Dcy.get4vnumu()
        P_numu[0] = P_numu[0]/1000.
        P_numu[1] = P_numu[1]/1000.
     
        P_numu4v  = T4V.TLorentzVector(P_numu[1][0], \
                                       P_numu[1][1], \
                                       P_numu[1][2], \
                                       P_numu[0] )

        if NeutrinoEventInstance.__Debug:
            print("            ----> Rest frame P_numu4v (GeV):",
                  P_numu4v.E(), P_numu4v.Px(), P_numu4v.Py(), P_numu4v.Pz())

        P_numu4v.Boost(b)

        if NeutrinoEventInstance.__Debug:
            print("            ----> nuSTORM frame P_numu4v (GeV):", 
                  P_numu4v.E(), P_numu4v.Px(), P_numu4v.Py(), P_numu4v.Pz())

        P_e[0]    = P_e4v.E()
        P_e[1][0] = P_e4v.Px()
        P_e[1][1] = P_e4v.Py()
        P_e[1][2] = P_e4v.Pz()
        
        P_nue[0]    = P_nue4v.E()
        P_nue[1][0] = P_nue4v.Px()
        P_nue[1][1] = P_nue4v.Py()
        P_nue[1][2] = P_nue4v.Pz()
        
        P_numu[0]    = P_numu4v.E()
        P_numu[1][0] = P_numu4v.Px()
        P_numu[1][1] = P_numu4v.Py()
        P_numu[1][2] = P_numu4v.Pz()
        
        return P_e, P_nue, P_numu
    

#--------  get/set methods:
    def getpmu(self):
        return deepcopy(self._pmu)

    def getTraceSpaceCoord(self):
        return deepcopy(self._TrcSpcCrd)

    def getpmuGen(self):
        return deepcopy(self._pmuGen)
    
    def getPb(self):
        Pb = np.array([self._pmuGen*self._pmuDirCos[0], \
                       self._pmuGen*self._pmuDirCos[1], \
                       self._pmuGen*self._pmuDirCos[2]  ])
        return Pb
    
    def gete4mmtm(self):
        return deepcopy(self._P_e)

    def getnue4mmtm(self):
        return deepcopy(self._P_nue)

    def getnumu4mmtm(self):
        return deepcopy(self._P_numu)
