#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class plane:
===========+

  Models a plane downstream from the nuStorm production straight

  Version: 2.2                                  Date:11 March 2021
  Author: Paul Kyberd
  
  Version 2.5   findHitPostion name clash - not resolved by different parameters.
                create two routines findHitPositionMuEvt and findHitPositionPiEvt
  
  Version 2.4   Add a findHitPosition for the pion decays

  Version 2.3   Production straight and flux plane from data file

  Version 2.2   Length of production straight past into constructor

  Version 2.1   remove the electron from consideration and do both neutrinos

  Version 2.0   Allow the x and y positions of the decaying muon to be off axis

  Version 1.0   Get the x, y position of the intersection of a neutrino trajectory
                with a plane at a distance sP (m) down stream of the end of the straight
                given the muonDecay instance


>>>>> this to modify 
  Class attributes:
  -----------------
  --np       : numpy class
      
  Instance attributes:
  --------------------
  zPos: postion of plane downstream of the straight
  particleSpecies: type of particle for which this plane was created 
    
  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : a plane instance
      __repr__ : One liner with call.
      __str__  : Dump of values of the plane

  Get/set methods:
    
  
  General methods:
    findHitPosition      : finds and stores the hit position of a particle with the plane

"""

from copy import deepcopy
import math
import numpy as np
import matplotlib.pyplot as plt
import sys
import nuSTORMPrdStrght as nuStrt

    
class plane:

    __Debug  = False


#--------  "Built-in methods":
    def __init__(self, straightLength, planePosition):

        self.zPos = planePosition + straightLength
        return

    def __repr__(self):
        return "plane(zPos, particleSpecies)\n"\
        "          z is distance downstream of the start of the nuStorm Straight\n" \
        "          Particle Species is corresponding particle for the plane "


    def __str__(self):
        return "plane : z(m) = %g \r\n" \
                % (self.zPos)

# Find the intercept of the particles with plane and store it along with the energy for muon events
    def findHitPositionMuEvt(self, nuEvt):
#  Decay point information
        DecayPntX = nuEvt.getTraceSpaceCoord()[1]
        DecayPntY = nuEvt.getTraceSpaceCoord()[2]
        DecayPntZ = nuEvt.getTraceSpaceCoord()[3]

#  Now the information for both neutrinos
        eE = nuEvt.getnue4mmtm()[0]         # number
        pE = nuEvt.getnue4mmtm()[1]         # 3 vector
        eMu = nuEvt.getnumu4mmtm()[0]
        pMu = nuEvt.getnumu4mmtm()[1]

#  Distance to the plane
        deltaZ = self.zPos - DecayPntZ
        if self.__Debug: print ("path length is ", deltaZ)
#  Position of the nue hit               # need to check by hand - no test

        hitE=[]
        xPnt = (DecayPntX + pE[0]*deltaZ/pE[2])
        yPnt = (DecayPntY + pE[1]*deltaZ/pE[2])
        hitE.append(xPnt)                                 # x
        hitE.append(yPnt)                                 # y
        hitE.append(self.zPos)                            # z
        hitE.append(np.sqrt(xPnt*xPnt + yPnt*yPnt))       # R
        hitE.append(math.atan2(yPnt, xPnt))               # phi
        hitE.append(pE[0])                               # px
        hitE.append(pE[1])                               # py
        hitE.append(pE[2])                               # pz
        hitE.append(eE)                                   # E

        if self.__Debug: print (" hit: x, y, R, phi, E ", hitE[0], " ", hitE[1], " ", hitE[2], " ", hitE[3], " ", hitE[4])

#  Position of the numu hit               # need to check by hand - no test
        hitMu=[]
        xPnt = (DecayPntX + pMu[0]*deltaZ/pMu[2])
        yPnt = (DecayPntY + pMu[1]*deltaZ/pMu[2])
        hitMu.append(xPnt)                                 # x
        hitMu.append(yPnt)                                 # y
        hitMu.append(self.zPos)                            # z
        hitMu.append(np.sqrt(xPnt*xPnt + yPnt*yPnt))       # R
        hitMu.append(math.atan2(yPnt, xPnt))               # phi
        hitMu.append(pMu[0])                               # px
        hitMu.append(pMu[1])                               # py
        hitMu.append(pMu[2])                               # pz
        hitMu.append(eMu)                                  # E


        if self.__Debug: print (" hit: x, y, R, phi, E ", hitMu[0], " ", hitMu[1], " ", hitMu[2], " ", hitMu[3], " ", hitMu[4])

        return hitE, hitMu

# Find the intercept of the particle with plane and store it along with the energy for pion flash events

    def findHitPositionPiEvt(self, piEvt):
#  Decay point information
        DecayPntX = piEvt.getTraceSpaceCoord()[1]
        DecayPntY = piEvt.getTraceSpaceCoord()[2]
        DecayPntZ = piEvt.getTraceSpaceCoord()[3]

#  Now the information for the neutrino
        eMu = piEvt.getnumu4mmtm()[0]       # number
        pMu = piEvt.getnumu4mmtm()[1]       # 3 vector

#  Distance to the plane
        deltaZ = self.zPos - DecayPntZ
#  Position of the numu hit               # need to check by hand - no test
        hitMu=[]
        xPnt = (DecayPntX + pMu[0]*deltaZ/pMu[2])
        yPnt = (DecayPntY + pMu[1]*deltaZ/pMu[2])
        if (xPnt < -200.0):
            print ("xPnt is ", xPnt, "   yPnt is ", yPnt)
            print ("piEvt is ", piEvt)
        hitMu.append(xPnt)                                 # x
        hitMu.append(yPnt)                                 # y
        hitMu.append(self.zPos)                            # z
        hitMu.append(np.sqrt(xPnt*xPnt + yPnt*yPnt))       # R
        hitMu.append(math.atan2(yPnt, xPnt))               # phi
        hitMu.append(pMu[0])                               # px
        hitMu.append(pMu[1])                               # py
        hitMu.append(pMu[2])                               # pz
        hitMu.append(eMu)                                  # E


        if self.__Debug: print (" hit: x, y, R, phi, E ", hitMu[0], " ", hitMu[1], " ", hitMu[2], " ", hitMu[3], " ", hitMu[4])

        return hitMu
    def findHitPositionPiFlash(self, nuFlash):

#  Decay point information
        DecayPntX = nuFlash.x()
        DecayPntY = nuFlash.y()
        DecayPntZ = nuFlash.z()
        pMu = nuFlash.p()[1]
        eMu = nuFlash.p()[0]

#  Distance to the plane
        deltaZ = self.zPos - DecayPntZ
        if self.__Debug: print ("path length is ", deltaZ)

#  Position of the numu hit               # need to check by hand - no test
        hitMu=[]
        xPnt = (DecayPntX + pMu[0]*deltaZ/pMu[2])
        yPnt = (DecayPntY + pMu[1]*deltaZ/pMu[2])
        if (xPnt < -200.0):
            print ("xPnt is ", xPnt, "   yPnt is ", yPnt)
            print ("nuFlash is ", nuFlash)
        hitMu.append(xPnt)                                 # x
        hitMu.append(yPnt)                                 # y
        hitMu.append(self.zPos)                            # z
        hitMu.append(np.sqrt(xPnt*xPnt + yPnt*yPnt))       # R
        hitMu.append(math.atan2(yPnt, xPnt))               # phi
        hitMu.append(pMu[0])                               # px
        hitMu.append(pMu[1])                               # py
        hitMu.append(pMu[2])                               # pz
        hitMu.append(eMu)                                  # E


        if self.__Debug: print (" hit: x, y, R, phi, E ", hitMu[0], " ", hitMu[1], " ", hitMu[2], " ", hitMu[3], " ", hitMu[4])

        return hitMu
