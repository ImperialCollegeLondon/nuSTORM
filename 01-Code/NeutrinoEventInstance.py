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
  _Detection_P_e    : Checking if the electron is detected 
  _Detection_P_nue  : Checking if the electron neutrino is detected
  _Detection_P_numu : Checking if the muon neutrino is detected
      
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
    is_e_detected     : Returns True if the electron is detected
    is_nue_detected   : Returns True if the electron neutrino is detected
    is_numu_detected  : Returns True if the muon neutrino is detected 

  General methods:
    CreateNeutrinos      : Manager for neutrino decay, returns z (m) of decay (real), P_e, P_nue, P_numu, RestFrame
                           Restframe contains a dump of the instance attributes of the MuonDecay class
    GenerateDcyPhaseSpace: Generates trace-space position of decay
    GenerateLongiPos     : Returns s, z of decay
    BeamDir              : Returns position and momentum of the beam and rotation operator correspoinding to the direction of velocity 
    Boost2nuSTORM        : Boost to nuSTORM rest frame -- i.e. boost to pmu
    RotMat               : Returns a rotation operator matrix that rotates a (0,0,z) vector along the direction of exact muon momentum
    RotnBoost            : Operator; rotates and boosts rest-frame coordinates to nuSTORM frame
    isDetected           : Returns True if the decay product is detected, otherwise False.

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
import math

muCnst = MuonConst.MuonConst()

global debug
global Debug
debug=0
Debug=0

    
class NeutrinoEventInstance:

    __mumass = muCnst.mass()/1000.
    __sol    = muCnst.SoL()

    __Debug  = False

#--------  "Built-in methods":
    def __init__(self, pmu=5., filename=None):

        nuStrt = nuPrdStrt.nuSTORMPrdStrght(filename)

        self._pmu = pmu
        self._ct, self._TrcSpcCrd, self._pmuGen, self._P_e, self._P_nue, self._P_numu, self._Pb, self._Detection_P_e, self._Detection_P_nue, self._Detection_P_numu = self.CreateNeutrinos(nuStrt)

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
        Circumference  = nuStrt.Circumference()
        ArcLen         = nuStrt.ArcLen()       
        ArcRad         = ArcLen / math.pi   #.. KL mod

        #debug=0            #.. Set debug = 1 to debug this method

        if(debug==1):
           print("--------Begining OF INSTANCE-------")

        #.. Prepare--get neutrino decay instance in muon rest frame:
        z = 2.* PrdStrghtLngth
        Dcy = 0
        if NeutrinoEventInstance.__Debug:
            print("NeutrinoEventInstance.CreateNeutrinos: find valid decay")
        while z > (PrdStrghtLngth+ArcRad+1.):    #.. KL mod
            if isinstance(Dcy, MuonDecay.MuonDecay):
                del Dcy
            
            Pmu0 = self.getpmu()
            Pmu  = nuStrt.GenerateMmtm(Pmu0)
            Emu   = np.sqrt(Pmu**2 + NeutrinoEventInstance.__mumass**2)
            beta  = Pmu / Emu
            gamma = Emu / NeutrinoEventInstance.__mumass
            v    = beta * NeutrinoEventInstance.__sol
           # Tmax = 10*nuStrt.Circumference() / (gamma * v)
            Dcy = MuonDecay.MuonDecay()                 # For no truncation of lifetime distribution
           # Dcy = MuonDecay.MuonDecay(Tmax=Tmax)         # Truncated lifetime distribution
            DcyCoord, pmuGen = self.GenerateDcyPhaseSpace(Dcy, Pmu, nuStrt)
            z  = DcyCoord[3]
            ct = Dcy.getLifetime()
            s  = DcyCoord[0]
            Pb = self.BeamDir(s, nuStrt, Pmu)[4]
            Pm = np.array([DcyCoord[4], DcyCoord[5], DcyCoord[6]])
           
            if(debug==1):
                  print("CreateNeutrinos: Pmu", Pmu)
                  print("CreateNeutinos: Pb", Pb)
                  print("CreateNeutinos: Pm", Pm)
       
        if z > (PrdStrghtLngth+ArcRad+1.):
            print("NeutrinoEvenInstance.CreateNeutrinos Alarm:", z)
        if NeutrinoEventInstance.__Debug:
            print("NeutrinoEventInstance.CreateNeutrinos: decay at z =", z)
            print("----> Dcy:", Dcy.__str__())
            
        #.. Boost to nuSTORM frame:
        if NeutrinoEventInstance.__Debug:
            print("NeutrinoEventInstance.CreateNeutrinos: rotate and boost to nuSTORM rest frame:")
        P_e, P_nue, P_numu = self.Boost2nuSTORM(Dcy, nuStrt,s, Pmu, Pm)
        if NeutrinoEventInstance.__Debug:
            print("----> P_e   :", P_e)
            print("----> P_nue :", P_nue)
            print("----> P_numu:", P_numu)

        R = [DcyCoord[1],DcyCoord[2],DcyCoord[3]]

        Detection_P_e    = self.isDetected(R, P_e[1], nuStrt)
        Detection_P_nue  = self.isDetected(R, P_nue[1], nuStrt)
        Detection_P_numu = self.isDetected(R, P_numu[1], nuStrt)


        del Dcy
        return ct, DcyCoord, pmuGen, P_e, P_nue, P_numu, Pb, Detection_P_e, Detection_P_nue, Detection_P_numu

#.. Trace space coordinate generation: array(s, x, y, z, x', y')
    def GenerateDcyPhaseSpace(self, Dcy, Pmu, nuStrt):
        #Debug=0  #.. Set Debug=1 to debug this method
        if(Debug==1):
              print('GenerateDcyPhaseSpace starting')
        coord = np.array([0., 0., 0., 0., 0., 0., 0.])

        #.. longitudinal position, "s", z:
        coord[0] = self.GenerateLongiPos(Dcy, Pmu)                                #s
       
        theta=self.BeamDir(coord[0], nuStrt, Pmu)[3]      
     
        x, y, xp, yp = nuStrt.GenerateTrans(coord[0])
        
        #Test: Amplify the value of transverse momentum so they can be seen in the momentum plot 
        #xp=500*xp
        
        coord[1] = x*math.cos(theta)+self.BeamDir(coord[0], nuStrt, Pmu)[2][0]    #xm
       
        coord[2] = y+self.BeamDir(coord[0], nuStrt, Pmu)[2][1]                    #ym 
  
        coord[3] = -x*math.sin(theta)+self.BeamDir(coord[0], nuStrt, Pmu)[2][2]   #zm
      
        coord[4] = self.BeamDir(coord[0], nuStrt, Pmu)[4][0]  + xp*math.cos(theta) #pmx 
        
        coord[5] = self.BeamDir(coord[0], nuStrt, Pmu)[4][1]  + yp               #pmy
        
        coord[6] = self.BeamDir(coord[0], nuStrt, Pmu)[4][2]  - xp*math.sin(theta)  #pmz
        
        return coord, Pmu

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

        numberofturns  = int(s/Circumference)
        #print("Number of turns",numberofturns)
        where          = s%Circumference
        r              = ArcLen/math.pi #r = ArcRad
        
        #Debug=0   # Set Debug=1 to debug this method

        if ( PrdStrghtLngth>=where):
            if(Debug==1):
             print("The muon is in the production straight")

            theta = 0     #Angle with respect to Z axis


            BeamPos = [0,0,where]      # Position of the beam when the decay occurs
            Pb = [Pmu* math.sin(theta), 0, Pmu* math.cos(theta)]  #Longitudinal momentum of the muon when the decay occurs
            #Beam Momentum
    
            if(Debug==1):
              print('BeamDir: theta,s,Pb', theta, s, Pb)       
 
        if (PrdStrghtLngth+ArcLen>=where>PrdStrghtLngth):
            if(Debug==1):
             print("The muon is in the first arc")

            ArcLenCovered = where - PrdStrghtLngth
            theta= math.pi*ArcLenCovered/ArcLen #in radians
            
            BeamPos=[(r-r* math.cos(theta)),0,PrdStrghtLngth+r*math.sin(theta)]
            
            Pb=[Pmu* math.sin(theta), 0, Pmu* math.cos(theta)]

            if(Debug==1): 
              print('BeamDir: theta,s,Pb', theta, s, Pb)
 
        if (2*PrdStrghtLngth+ArcLen>=where>PrdStrghtLngth+ArcLen):
            if(Debug==1):
             print("The muon is in the return straight")

            theta=math.pi

            BeamPos=[2*r,0,PrdStrghtLngth-(where-ArcLen-PrdStrghtLngth)]
            Pb=[Pmu* math.sin(theta), 0, Pmu* math.cos(theta)]

            if(Debug==1):
              print('BeamDir: theta,s,Pb', theta, s, Pb)

        if (2*PrdStrghtLngth+2*ArcLen>=where>=2*PrdStrghtLngth+ArcLen):
            if(Debug==1):
             print("The muon is in the second/return arc")

            ArcLenCovered=where - 2*PrdStrghtLngth - ArcLen
            theta= math.pi + math.pi*ArcLenCovered/ArcLen #angle with respect to  Z axis
            
            BeamPos=[(r-r*math.cos(theta)),0,r*math.sin(theta)]
            Pb=[Pmu* math.sin(theta), 0, Pmu* math.cos(theta)]

            if(Debug==1):
              print('BeamDir: theta, s, Pb', theta, s, Pb)
        
        # dummy variables to preserve the ordering of the outputs 
        Dummy1=1
        Dummy2=1
          
        return Dummy1, Dummy2, BeamPos, theta, Pb


#.. Boost from muon rest frame to nuSTORM frame:
    def Boost2nuSTORM(self, Dcy, nuStrt, s, Pmu, Pm):
        ''' Present approximation is muon propagates along z axis, so, boost only 
            in preparation for later, include rotation matrix to transform from
            nustorm frame to frame with z axis along muon momentum and back '''
        #Debug=0  #.. Set Debug=1 to debug this method
       
        if(Debug==1):
              print('Boost2nuSTORM starting')
              
        Emu = np.sqrt(Pmu**2 + \
                                              NeutrinoEventInstance.__mumass**2)
        beta   = Pmu / Emu
        gamma  = Emu / NeutrinoEventInstance.__mumass
       
        theta= self.BeamDir(s, nuStrt, Pmu)[3]##!!!!!  s
        
        R = self.RotMat(Pm)[0]
        Rinv = self.RotMat(Pm)[1]

        if(Debug==1):
              print('Boost2nuSTORM: theta,s', theta, s)
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
     
        P_e    = self.RotnBoost(P_e, R, Rinv, theta, gamma, beta)
        if NeutrinoEventInstance.__Debug:
            print("NeutrinoEventInstance.Boost2nuSTORM: nuSTORM frame P_e:")
            print("----> P_e (GeV):", P_e)
            print("----> Dcy.get4ve:", Dcy.get4ve())
        
        P_nue    = Dcy.get4vnue()
        P_nue[0] = P_nue[0]/1000.
        P_nue[1] = P_nue[1]/1000.
    
        P_nue    = self.RotnBoost(P_nue, R, Rinv, theta,  gamma, beta)
        
        P_numu    = Dcy.get4vnumu()
        P_numu[0] = P_numu[0]/1000.
        P_numu[1] = P_numu[1]/1000.
     
        P_numu    = self.RotnBoost(P_numu, R, Rinv, theta, gamma, beta)
        
        return P_e, P_nue, P_numu
        
    def RotnBoost(self, P, R, Rinv, theta, gamma, beta):
        if NeutrinoEventInstance.__Debug:
            print("NeutrinoEventInstance.RotnBoos:")
            print("----> P, R, Rinv, gamma, beta:", P, "\n", \
                  R, "\n", Rinv, "\n", gamma, beta)
        #Debug=0  #.. Set Debug=1 to debug this method
        
        if(Debug==1):
              print('RotnBoost starting')
              print('Gamma; Beta', gamma, beta)
              print("P[1], 3 vector before boost", P[1])
            
      #  p3 = np.dot(Rinv, P[1])                  #Not necessary in this approximation
        #print("p3, 3 vector after rotation", p3) 
        
        Ec = P[0]
        Pc = P[1][2]
        
        Ef = gamma * (Ec + beta * Pc)
        Pf = gamma * (Pc + beta * Ec)

        p3 = []
       
        p3.append(P[1][0])
       
        p3.append( P[1][1])
        p3.append(Pf)
        p3 = np.array(p3)
       
        if(Debug==1):
          print("RotnBoost:p3, 3 vector after boost", p3)

        p3=np.dot(R,p3)

        if(Debug==1):
          print("RotnBoost:p3, 3 vector after rotation", p3)
        if(Debug==1):
          print("RotnBoost: theta", theta)
              
        Po    = [0., np.array([0., 0., 0.])]
        Po[0] = Ef
        Po[1] = p3

        if(Debug==1):
          print('Po',Po)

        return Po

##..Checking if a decay product is detected        Needs 
    def isDetected(self, Rm, P, nuStrt):

        HallWallDist   = nuStrt.HallWallDist()        
        DetHlfWdth     = nuStrt.DetHlfWdth()      
        DetLngth       = nuStrt.DetLngth()
        Hall2Det       = nuStrt.Hall2Det()
        PrdStrghtLngth = nuStrt.ProdStrghtLen()
 
        xdmax        =  DetHlfWdth
        xdmin        = -DetHlfWdth
        ydmax        =  DetHlfWdth                                  
        ydmin        = -DetHlfWdth
        zd           =  PrdStrghtLngth + HallWallDist + Hall2Det
        
        x0 = Rm[0]
        y0 = Rm[1]  
        z0 = Rm[2]

        l  = P[0]
        m  = P[1]  
        n  = P[2]

        t  =(zd - z0)/n   
        xd = x0 + t*l      
        yd = y0 + t*m
        
        if (t>0):
          if (xdmax>=xd>=xdmin) and (ydmax>=yd>=ydmin):
            return True
        else:                  
            return False

#.. Rotaion operator that rotates a (0,0,z) vector to exact muon direction (px,py,pz)
    def RotMat(self, P):
        Phat = P/np.linalg.norm(P)
        zhat = np.array([0,0,1])
        phi = math.acos(np.dot(zhat,Phat))
        n  = np.cross(zhat,Phat)
        if np.all(n==0):
                nhat = n
        else:
                nhat = n/np.linalg.norm(n)
        nx = nhat[0]
        ny = nhat[1]
        nz = nhat[2]
        R   =  np.array([[math.cos(phi) + nx**2*(1-math.cos(phi)),        #Rxx
                          nx*ny*(1-math.cos(phi)) - nz* math.sin(phi),    #Rxy 
                          nx*nz*(1-math.cos(phi)) + ny* math.sin(phi)],   #Rxz
                         [nx*ny*(1-math.cos(phi)) + nz* math.sin(phi),    #Ryx
                          math.cos(phi) + ny**2*(1-math.cos(phi)),        #Ryy
                          ny*nz*(1-math.cos(phi)) - nx* math.sin(phi)],   #Ryz
                         [nx*nz*(1-math.cos(phi)) - ny* math.sin(phi),    #Rzx
                          ny*nz*(1-math.cos(phi)) + nx* math.sin(phi),    #Rzy
                          math.cos(phi) + nz**2*(1-math.cos(phi))]])      #Rzz 
        Rinv = np.linalg.inv(R)
        return R, Rinv
 
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
    
    def getPb(self):
        return deepcopy(self._Pb)
    
    def is_e_detected(self):
        return deepcopy(self._Detection_P_e)
    
    def is_nue_detected(self):
        return deepcopy(self._Detection_P_nue)
    
    def is_numu_detected(self):
        return deepcopy(self._Detection_P_numu)
