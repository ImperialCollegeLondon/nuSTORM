#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model for calculating normalised numbers
========================================

    Assumes that nuSim code is in python path.

    @file eventRate.py

    @brief   calculates the event rate in the detector

    @author  Paul Kyberd


    @version    2.0
    @date       05 June 2023
    @author     Paul Kyberd

    Based on mkNuSpectra

    change the initialisation of plane to allow the definition of the position in 
    x,y and z
    @version    1.4
    @date       27 July 2022
    @author     Paul Kyberd

    Add differentiation between pion decay phase space and pion phase space at target;
    Add right muon beam momentum to event history;
    Add central muon momentum to introduce muon momentum acceptance cut;
    Add argument parser for easier handling when running simulation as batch jobs;
    Add correct momentum calculation from angles and absolute momentum;
    Add nuSTORM constant handling through nuSTORM constant class;
    @version    1.3
    @date       09 June 2022
    @author     Marvin Pfaff

    Add the event history at the end of the production straight
    @version     1.2
    @date        07 January 2022


    Add python logging
    @version     1.1
    @date        07 January 2022


    @version     1.0
    @date        05 October 2021

"""
#
# Class to do the calculation of the event rate normalisation
#
import os, sys
from datetime import datetime
from copy import deepcopy
import argparse
import numpy as np
import math as math
import logging
import PionConst as PC
import MuonConst as MC
import nuSTORMConst
import control
import histoManager
import nuSTORMPrdStrght as nuPrdStrt
import nuSTORMTrfLineCmplx as nuTrfLineCmplx
import PionEventInstance as piEvtInst
import NeutrinoEventInstance as nuEvtInst
import RandomGenerator as Rndm
import plane as plane
import particle as particle
import eventHistory as eventHistory

class normalisation:

    __version__ = 1.3

    def __init__(self, hFlag, muonMom=0):
        self._tlDcyCount = 0
        self._byndPSCount = 0
        self._PSDcyCount = 0
        self._muDcyCount = 0
        self._tlAngle = tlCmplxAngle*math.pi/180.0
        self._sth = math.sin(self._tlAngle)
        self._cth = math.cos(self._tlAngle)
        self._mup0 = muonMom
        self.__history = hFlag

# transform x,y,z and px,py,pz co-ordinates from the transfer line local co-ordinates to the
# global ones

    def tltoGlbl(self, xl, yl, zl, pxl, pyl, pzl):

        xg = xl - zl*self._sth
        yg = yl
        zg = zl*self._cth

        pxg = pxl*self._cth + pzl*self._sth
        pyg = pyl
        pzg = pzl*self._cth - pxl*self._sth

        return xg, yg, zg, pxg, pyg, pzg

# Deal with decay in the transfer line
    def tlDecay(self):

      if (self._tlDcyCount < printLimit):  
        print (" +++++ starting tlDecay +++++")
      self._tlDcyCount = self._tlDcyCount + 1
# decays at present just set the productionStraight particle with a weight of zero - but updated s and z (small value of pz
# for tracespace calculation)
      noParticle = particle.particle(runNumber, event, tlCmplxLength, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, "none")
      if (self.__history): eH.addParticle("productionStraight", noParticle)
# add a pion decay particle - set the time to the decay lifetime
      if (self._tlDcyCount < printLimit): print ("pi in tlDecay ", pi)
      dcytsc = pi.getTraceSpaceCoord()
      sd = dcytsc[0]
      xdl = dcytsc[1]
      ydl = dcytsc[2]
      zdl = dcytsc[3]
      xpd = dcytsc[4]
      ypd = dcytsc[5]
      pPion = pi.getppiGen()
      pzdl = np.sqrt(pPion**2/(1+xpd**2+ypd**2))
      pxdl = pzdl*xpd
      pydl = pzdl*ypd
      pi.getLifetime()
      td = lifetime*1E9 + t
      zdl = zdl - tlCmplxLength
# Now we need to move to the global co-ordinates from the transfer line local
      xd, yd, zd, pxd, pyd, pzd = self.tltoGlbl(xdl, ydl, zdl, pxdl, pydl, pzdl)
      pionTLDecay = particle.particle(runNumber, event, sd, xd, yd, zd, pxd, pyd, pzd, td, eventWeight, "pi+")
      if (self._tlDcyCount < printLimit): print ("tlDecay: about to add a particle")
      if (self.__history): eH.addParticle("pionDecay", pionTLDecay)
      if (self._tlDcyCount < printLimit): print (" pion at pionDecay: TL ", pionTLDecay)
# add the pion flash neutrino
      numu = pi.getnumu4mmtm()
      pxnul = numu[1][0]
      pynul = numu[1][1]
      pznul = numu[1][2]
      xd, yd, zd, pxnu, pynu, pznu = self.tltoGlbl(xdl, ydl, zdl, pxnul, pynul, pznul)
      nuFlashTL = particle.particle(runNumber, event, sd, xd, yd, zd, pxnu, pynu, pznu, td, eventWeight, "numu")
      if (self.__history): eH.addParticle("piFlashNu", nuFlashTL)
      if (self._tlDcyCount < printLimit): print ("at piFlashNu: TL")
# add the muon from the pion flash ... not going to track this, weight non zero, but don't track further ?
      mu = pi.getmu4mmtm()
      pxmul = mu[1][0]
      pymul = mu[1][1]
      pzmul = mu[1][2]
      xd, yd, zd, pxmu, pymu, pzmu = self.tltoGlbl(xdl, ydl, zdl, pxmul, pymul, pzmul)
      muonProdTL = particle.particle(runNumber, event, sd, xd, yd, zd, pxmu, pymu, pzmu, td, eventWeight, "mu+")
      if (self.__history): eH.addParticle("muonProduction",muonProdTL)
      if (self._tlDcyCount < printLimit): print ("at muonProduction: TL")
# extrapolate the neutrino the the detector plane
      hitMu = fluxPlane.findHitPositionPiFlash(nuFlashTL)
      if (self._tlDcyCount < printLimit): print ("hit position of neutrino ", hitMu)
#  fill the event History
      numuX = hitMu[0] + detectorPosition[0]
      numuY = hitMu[1] + detectorPosition[1]
      numuZ = hitMu[2] + detectorPosition[2]
      dsNumu = math.sqrt((xd-numuX)**2 + (yd-numuY)**2 + (zd-numuZ)**2)
      sNumu = sd + dsNumu
      tNumu = td + dsNumu*1E9/c
      if (self._tlDcyCount < printLimit): print ( "sNumu is ", sNumu, "    dsNumu is ", dsNumu, "     tNumu is ", tNumu, "    c is ", c)
      if ((abs(numuX) < 100) and (abs(numuY) < 100)):
          eW = eventWeight
      else:
          eW = 0.0
      numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, pxnu, pynu, pznu, tNumu, eW, "numu")
      if (self.__history): eH.addParticle("numuDetector", numuDetector)
      if (self._tlDcyCount < printLimit): print ("numu at detector")

      return
#
# Deal with decay beyond the production straight -------------------------------------------------
#
    def beyondPS(self):

      if (self._byndDcyCount < printLimit):  
        print (" +++++ starting beyondPS +++++")
      self._byndPSCount = self._byndPSCount + 1
      # extraoplate the muon to the end of the production straight
      tscTrgtLcl = pi.getLclTraceSpaceCoord()
      xtl  = tscTrgtLcl[1]
      ytl  = tscTrgtLcl[2]
      xptl = tscTrgtLcl[4]
      yptl = tscTrgtLcl[5]
      pztl = np.sqrt(pPion**2/(1+xptl**2+yptl**2))
      pxtl = xptl*pztl
      pytl = yptl*pztl
      ttl  = t
      sEnd = tlCmplxLength + psLength
      zEnd = psLength
      xEnd = xtl #+ dZ*pxMu/pzMu
      yEnd = ytl #+ dZ*pyMu/pzMu
      dZ   = sEnd
      dFlown = math.sqrt((xEnd-xtl)*(xEnd-xtl) + (yEnd-ytl)*(yEnd-ytl) + dZ*dZ)
      piVel  = math.sqrt(pxtl*pxtl + pytl*pytl + pztl*pztl)*c/Epion
      tFlown = dFlown/piVel
      tEnd = ttl + tFlown
      piEnd = particle.particle(runNumber, event, sEnd, xEnd, yEnd, zEnd, pxtl, pytl, pztl, tEnd, eventWeight, "pi+")
      if (self.__history): eH.addParticle("prodStraightEnd", piEnd)
      #  add a pion decay particle - set the time to the decay lifetime and the s to the pathlength, eventweight to full
      # x,y,z to 0.0 and px,py to 0.0, pz to 0.01 so constructor does not
      dcytsc = pi.getTraceSpaceCoord()
      sd = dcytsc[0]
      xd = 0.0
      yd = 0.0
      zd = 0.0
      pxd = 0.0
      pyd = 0.0
      pzd = 0.01
      td = pi.getLifetime()*1E9 + t
      if (self._byndPSCount < printLimit): print ("pi in beyondPS: decayLength ", sd)
      pionLostDecay = particle.particle(runNumber, event, sd, xd, yd, zd, pxd, pyd, pzd, td, eventWeight, "none")
      if (self.__history): eH.addParticle("pionDecay", pionLostDecay)
      if (self._byndPSCount < printLimit):  print ("at pionDecay lost")
# add the pion flash neutrino ... set everything to zero - including eventWeight
      nuFlashLost = particle.particle(runNumber, event, sd, xd, yd, zd, pxd, pyd, pzd, td, 0.0, "none")
      if (self._byndPSCount < printLimit):  print ("at piFlashNu lost")
      if (self.__history): eH.addParticle("piFlashNu", nuFlashLost)
# add the muon from the pion flash ... set everything to zero - including eventWeight
      muonProdLost = particle.particle(runNumber, event, sd, xd, yd, zd, pxd, pyd, pzd, td, 0.0, "none")
      if (self._byndPSCount < printLimit):  print ("at muonProduction lost")
      if (self.__history): eH.addParticle("muonProduction",muonProdLost)
      return

##
# calculate where in the ring, the decay occurred given the path length in the ring
#
    def ring(self, pathLength):

      arcRadius = 50.0
      sectionLengths=[180.0, arcRadius*math.pi, 180.0, arcRadius*math.pi]

      sectionPnt = 0
      while pathLength > 0.0:
        pathLength = pathLength - sectionLengths[sectionPnt]
        sectionPnt = (sectionPnt+1)%4

      sectionPnt = (sectionPnt-1)%4
      decayPnt = pathLength + sectionLengths[sectionPnt]
      if (sectionPnt == 0):
        theta = 0.0
        zpos = 0.0 + decayPnt
        deltaX = 0.0
        direction = 0.0
      elif (sectionPnt == 2):
        theta = -math.pi
        zpos = 180.0 - decayPnt
        deltaX = -2*arcRadius
        direction = math.pi
      elif (sectionPnt == 1):
        theta = decayPnt/arcRadius - 0.5*math.pi
        zpos = 180.0 + arcRadius*math.cos(theta)
        deltaX = 0.0 - arcRadius - arcRadius*math.sin(theta)
        direction = 2.0*math.pi - decayPnt/arcRadius
      elif (sectionPnt == 3):
        theta = decayPnt/arcRadius + 0.5*math.pi
        zpos = 0.0 + arcRadius*math.cos(theta)
        deltaX = -arcRadius - arcRadius*math.sin(theta)
        direction = math.pi - decayPnt/arcRadius

      return deltaX, zpos, direction

#
# Deal with decay in the production straight -------------------------------------------------
#
    def decayPiInPS(self):

      if (self._PSDcyCount < printLimit):  
        print (" +++++ starting decayPiInPS +++++")
      self._PSDcyCount = self._PSDcyCount + 1
      dcytsc = pi.getTraceSpaceCoord()
      sd = dcytsc[0]
      xd = dcytsc[1]
      yd = dcytsc[2]
      xpd = dcytsc[4]
      ypd = dcytsc[5]
#  need the lifetime in the nuStorm frame in ns
      piLifetime = pi.getLifetime()*1E9*Epion/(piMass)
      td = piLifetime + t
      zd = sd - tlCmplxLength
      pzd = np.sqrt(pPion**2/(1+xpd**2+ypd**2))
      pxd = pzd*xpd
      pyd = pzd*ypd
      pionPSDecay = particle.particle(runNumber, event, sd, xd, yd, zd, pxd, pyd, pzd, td, eventWeight, "pi+")
      if (td < 150.0):
        print ("error in the time")
#      hLifetime.Fill(piLifetime)

      if (self.__history): eH.addParticle("pionDecay", pionPSDecay)
      if (self._PSDcyCount < printLimit): print ("pionDecay in PS")
# add the muon
      mu = pi.getmu4mmtm()
      pxMu = mu[1][0]
      pyMu = mu[1][1]
      pzMu = mu[1][2]
      eMu = mu[0]
      muProd = particle.particle(runNumber, event, sd, xd, yd, zd, pxMu, pyMu, pzMu, td, eventWeight, "mu+")
      if (self.__history): eH.addParticle("muonProduction", muProd)
# extraoplate the muon to the end of the production straight
      dZ = psLength - zd
      zEnd = psLength
      xEnd = xd #+ dZ*pxMu/pzMu
      yEnd = yd #+ dZ*pyMu/pzMu
      sEnd = tlCmplxLength + psLength
      dFlown = math.sqrt((xEnd-xd)*(xEnd-xd) + (yEnd-yd)*(yEnd-yd) + dZ*dZ)
      muVel = math.sqrt(pxMu*pxMu + pyMu*pyMu + pzMu*pzMu)*c/eMu
      tFlown = dFlown/muVel
      tEnd = td + tFlown
      muTSC = nuEvt.getTraceSpaceCoord()
      sDcy = muTSC[0]
      if (sDcy > tlCmplxLength+psLength):
          muEnd = particle.particle(runNumber, event, sEnd, xEnd, yEnd, zEnd, pxMu, pyMu, pzMu, tEnd, eventWeight, "mu+")
          if (self.__history): eH.addParticle("prodStraightEnd", muEnd)
      else:
          noParticle = particle.particle(runNumber, event, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0,   "none")
          if (self.__history): eH.addParticle("prodStraightEnd", noParticle)

      if (self._PSDcyCount < printLimit): print ("muonProduction in production straight")
# add the pion flash neutrino
      numu = pi.getnumu4mmtm()
      pxnu = numu[1][0]
      pynu = numu[1][1]
      pznu = numu[1][2]
      eNu = numu[0]

      nuFlash = particle.particle(runNumber, event, sd, xd, yd, zd, pxnu, pynu, pznu, td, eventWeight, "numu")
      if (self.__history): eH.addParticle("piFlashNu", nuFlash)
      if (self._PSDcyCount < printLimit):  
        print ("piFlashNu PS")
        print ("nuFlash is ", nuFlash)
# extrapolate the neutrino the the detector plane
      if FlshAtDetFlg:
          hitMu = fluxPlane.findHitPositionPiFlash(nuFlash)
          if (self._PSDcyCount < printLimit): print ("hit position of neutrino ", hitMu)
#  fill the event History
          numuX = hitMu[0] + detectorPosition[0]
          numuY = hitMu[1] + detectorPosition[1]
          numuZ = hitMu[2] + detectorPosition[2]
#   Because the hit is in detector local we need to transform the machine local for the calculation.
#   At present the x and y are the same but the z is at 230
          dsNumu = math.sqrt((numuX-xd)**2 + (numuY-yd)**2 + ((numuZ+230.0)-zd)**2)
          sNumu = sd + dsNumu
          tNumu = td + dsNumu*1E9/c
          if (self._PSDcyCount < printLimit): print ( "xdiff is  ", xd-numuX, "    ydiff is ", yd-numuY, "     zdiff is ", (numuZ+230.0)-zd, "   numuX is ", numuX, "    numuY is ", numuY, "   numuZ is ", numuZ)
          if (self._PSDcyCount < printLimit): print ( "sNumu is ", sNumu, "    dsNumu is ", dsNumu, "     tNumu is ", tNumu, "    c is ", c)
          if ((abs(numuX) < 10.0) and (abs(numuY) < 10.0)):
            eW = eventWeight
            ENumuSpectra = math.sqrt(pxnu*pxnu + pynu*pynu + pznu*pznu)
            hNumuDet.Fill(ENumuSpectra)        
            numuDetTime.Fill(tNumu)  
          else:
            eW = 0.0
          numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, pxnu, pynu, pznu, tNumu, eW, "numu")
          if (self.__history): eH.addParticle("numuDetector", numuDetector)
          if (self._PSDcyCount < printLimit): print ("numu at detector", numuDetector)

#
# Muon decays -------------------------------------------------
#
    def decayMuons(self):

# ... We can allow any muon which is created in the production straight to decay without worrying about
# ... acceptance

      if (self._PSDcyCount < printLimit):  
        print (" +++++ starting decayMuons")
      piTraceSpaceCoord = pi.getTraceSpaceCoord()
      mucostheta = pi.getcostheta()
      mu4mom = pi.getmu4mmtm()
      if (self._muDcyCount < printLimit):
            print ("piTraceSpaceCoord is ", piTraceSpaceCoord)
            print ("muon 4 momentum is ", mu4mom)
            print ('muoncostheta is ' , mucostheta)
      mup0 = self._mup0
      if self._mup0 == 0.:
          print("ERROR: Muon central momentum is 0.! Interrupting simulation...")
          exit()
      Absorbed = nuEvt.Absorption(piTraceSpaceCoord, mu4mom, mucostheta, mup0)
#        if (Absorbed == False): print ("absorbed is ", Absorbed, "  for event ", event)
#  if the muon doesn't make it in the ring acceptance ... it is absorbed ... so for the
#  muon we create a suitable muon - but the other particles are all put to null values so
#  all values bar run and event are put to zero
      if (Absorbed):
        testParticle = particle.particle(runNumber, event, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0,   "none")
        if (self.__history): eH.addParticle("muonDecay", testParticle)
        if (self._muDcyCount < printLimit): print ("absorbed")
        if (PSMuonsFlag):
          muTSC = nuEvt.getTraceSpaceCoord()
          sDcy = muTSC[0]
#TL decay
#          if (sDcy < tlCmplxLength):
# ps decay
          if (sDcy < tlCmplxLength+psLength):
            print (f"========= muon decay in the production straight ============ {sDcy}")
# Muon Decay
            xDcy = muTSC[1]
            yDcy = muTSC[2]
            zDcy = muTSC[3]
            pxDcy = nuEvt.getPb()[0]
            pyDcy = nuEvt.getPb()[1]
            pzDcy = nuEvt.getPb()[2]
#   need to transform lifetimes to the ring frame
            ppi = pi.getppi()
            pmuon = pi.getmu4mmtm()
            piMass = piCnst.mass()/1000.0
            muMass = muCnst.mass()/1000.0
            gammaPi = math.sqrt(1 + (ppi*ppi)/(piMass*piMass))
            gammaMu = math.sqrt(1 + (pmuon[0]*pmuon[0])/(muMass*muMass))
            tDcy = (gammaPi*pi.getLifetime()+gammaMu*nuEvt.getLifeTime())*1E9 + t
            if (self._muDcyCount < printLimit): print ("decayMuons: ", "Pi lifetime is ", pi.getLifetime()*1E9, "   mu lifetime is ", nuEvt.getLifeTime()*1E9, "   t is ", t, "   tDcy is ", tDcy)
            muDecay = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, pxDcy, pyDcy, pzDcy, tDcy, eventWeight, "mu+")
            if (self.__history): eH.addParticle("muonDecay", muDecay)
            self._muDcyCount = self._muDcyCount + 1
            if (self._muDcyCount < printLimit): print ("muDecay is ", muDecay)

# electron production
            e4mmtm = nuEvt.gete4mmtm()
            e3mmtm = e4mmtm[1]
            eProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, e3mmtm[0], e3mmtm[1], e3mmtm[2], tDcy, eventWeight, 'e+')
            if (self.__history): eH.addParticle('eProduction', eProd)
            if (self._muDcyCount < printLimit): print ("eProduction is ", eProd)

# numu production
            numu4mmtm = nuEvt.getnumu4mmtm()
            numu3mmtm = numu4mmtm[1]
            numuPx = numu3mmtm[0]
            numuPy = numu3mmtm[1]
            numuPz = numu3mmtm[2]
            numuProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, numu3mmtm[0], numu3mmtm[1], numu3mmtm[2], tDcy, eventWeight, 'numuBar')
            if (self.__history): eH.addParticle('numuProduction', numuProd)
            if (self._muDcyCount < printLimit): print ("numuProduction is ", numuProd)

# nue production
            nue4mmtm = nuEvt.getnue4mmtm()
            nue3mmtm = nue4mmtm[1]
            nuePx = nue3mmtm[0]
            nuePy = nue3mmtm[1]
            nuePz = nue3mmtm[2]
            nueProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, nue3mmtm[0], nue3mmtm[1], nue3mmtm[2], tDcy, eventWeight, 'nue')
            if (self.__history): eH.addParticle('nueProduction', nueProd)
            if (self._muDcyCount < printLimit): print ("nueProduction is ", nueProd)
# and finally extrapolate to the neutrino detector
#   hitx is x, y, z, R, phi, px, py, pz, E
            hitE,hitMu=fluxPlane.findHitPositionMuEvt(nuEvt)
#  numu extrapolation to the detector - returned in detector co-ordinates so z = 0 - transform to rinhg co-ordinates 
            numuX = hitMu[0] + detectorPosition[0]
            numuY = hitMu[1] + detectorPosition[1]
            numuZ = hitMu[2] + detectorPosition[2]
            dsNumu = math.sqrt((xDcy-numuX)**2 + (yDcy-numuY)**2 + (zDcy-numuZ)**2)
            sNumu = sDcy + dsNumu
            tNumu = tDcy + dsNumu*1E9/c
            if (self._muDcyCount < printLimit): print ( "sNumu is ", sNumu, "    dsNumu is ", dsNumu, "     tNumu is ", tNumu, "     tDecay is ", tDcy)
            if ((abs(numuX) < 2.50) and (abs(numuY) < 2.50)):
                eW = eventWeight
                ENumuSpectra = math.sqrt(numuPx*numuPx + numuPy*numuPy + numuPz*numuPz)
                hNumuDet.Fill(ENumuSpectra)
            else:
                eW = 0.0
            numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, numuPx, numuPy, numuPz, tNumu, eW, "numu")
            if (self.__history): eH.addParticle("numuDetector", numuDetector)
            if (self._muDcyCount < printLimit): print ("numu at detector")
#  nue extrapolation to the detector
            nueX = hitE[0] + detectorPosition[0]
            nueY = hitE[1] + detectorPosition[1]
            nueZ = hitE[2] + detectorPosition[2]
            dsNue = math.sqrt((xDcy-nueX)**2 + (yDcy-nueY)**2 + (zDcy-nueZ)**2)
            sNue = sDcy + dsNue
            tNue = tDcy + dsNue*1E9/c
            if ((abs(nueX) < 2.50) and (abs(nueY) < 2.50)):
                eW = eventWeight
                ENueSpectra = math.sqrt(nuePx*nuePx + nuePy*nuePy + nuePz*nuePz)
                hNueDet.Fill(ENueSpectra)
            else:
                eW = 0.0
            nueDetector = particle.particle(runNumber, event, sNue, nueX, nueY, nueZ, nuePx, nuePy, nuePz, tNue, eW, "nue")
            if (self.__history): eH.addParticle("nueDetector", nueDetector)
            if (self._muDcyCount < printLimit): print ("nue at detector")

# finished dealing with muon decays in the production straight
      else:
        if (ringMuonsFlag):
            self._muDcyCount = self._muDcyCount + 1
            if (self._muDcyCount < printLimit): print ("========= muon decay in the ring ===========")
            muTSC = nuEvt.getTraceSpaceCoord()
            if (self._muDcyCount < printLimit): print ("muTSC is ", muTSC)
# Muon Decay
            sDcy = muTSC[0]
            xDcy = muTSC[1]
            yDcy = muTSC[2]
            zDcy = muTSC[3]
            pxDcy = nuEvt.getPb()[0]
            pyDcy = nuEvt.getPb()[1]
            pzDcy = nuEvt.getPb()[2]
#   need to transform lifetimes to the ring frame
            ppi = pi.getppi()
            pmuon = pi.getmu4mmtm()
            piMass = piCnst.mass()/1000.0
            muMass = muCnst.mass()/1000.0
            gammaPi = math.sqrt(1 + (ppi*ppi)/(piMass*piMass))
            gammaMu = math.sqrt(1 + (pmuon[0]*pmuon[0])/(muMass*muMass))
            tDcy = (gammaPi*pi.getLifetime()+gammaMu*nuEvt.getLifeTime())*1E9 + t
            if (self._muDcyCount < printLimit): print ("decayMuons: ", "Pi lifetime is ", pi.getLifetime()*1E9, "   mu lifetime is ", nuEvt.getLifeTime()*1E9, "   t is ", t)            
            muDecay = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, pxDcy, pyDcy, pzDcy, tDcy, eventWeight, "mu+")
            if (self.__history): eH.addParticle("muonDecay", muDecay)
            if (self._muDcyCount < printLimit): print ("muDecay is ", muDecay)

# electron production
            e4mmtm = nuEvt.gete4mmtm()
            e3mmtm = e4mmtm[1]
            eProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, e3mmtm[0], e3mmtm[1], e3mmtm[2], tDcy, eventWeight, 'e+')
            if (self.__history): eH.addParticle('eProduction', eProd)
            if (self._muDcyCount < printLimit): print ("eProduction is ", eProd)

# numu production
            numu4mmtm = nuEvt.getnumu4mmtm()
            numu3mmtm = numu4mmtm[1]
            numuPx = numu3mmtm[0]
            numuPy = numu3mmtm[1]
            numuPz = numu3mmtm[2]
            numuProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, numu3mmtm[0], numu3mmtm[1], numu3mmtm[2], tDcy, eventWeight, 'numuBar')
            if (self.__history): eH.addParticle('numuProduction', numuProd)
            if (self._muDcyCount < printLimit): print ("numuProduction is ", numuProd)

# nue production
            nue4mmtm = nuEvt.getnue4mmtm()
            nue3mmtm = nue4mmtm[1]
            nuePx = nue3mmtm[0]
            nuePy = nue3mmtm[1]
            nuePz = nue3mmtm[2]
            nueProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, nue3mmtm[0], nue3mmtm[1], nue3mmtm[2], tDcy, eventWeight, 'nue')
            if (self.__history): eH.addParticle('nueProduction', nueProd)
            if (self._muDcyCount < printLimit): print ("nueProduction is ", nueProd)

# and finally extrapolate to the neutrino detector
#   hitx is x, y, z, R, phi, px, py, pz, E
            hitE,hitMu=fluxPlane.findHitPositionMuEvt(nuEvt)
#  numu extrapolation to the detector
            numuX = hitMu[0] + detectorPosition[0]
            numuY = hitMu[1] + detectorPosition[1]
            numuZ = hitMu[2] + detectorPosition[2]
            dsNumu = math.sqrt((xDcy-numuX)**2 + (yDcy-numuY)**2 + (zDcy-numuZ)**2)
            sNumu = sDcy + dsNumu
            tNumu = tDcy + dsNumu*1E9/c
            if (self._muDcyCount < printLimit): print ( "sNumu is ", sNumu, "    dsNumu is ", dsNumu, "     tNumu is ", tNumu)
            if ((abs(numuX) < 2.50) and (abs(numuY) < 2.50)):
                eW = eventWeight
                ENumuSpectra = math.sqrt(numuPx*numuPx + numuPy*numuPy + numuPz*numuPz)
                hNumuDet.Fill(ENumuSpectra)
            else:
                eW = 0.0
            numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, numuPx, numuPy, numuPz, tNumu, eW, "numu")
            if (self.__history): eH.addParticle("numuDetector", numuDetector)
            if (self._muDcyCount < printLimit): print ("numu at detector")

#  nue extrapolation to the detector
            nueX = hitE[0] + detectorPosition[0]
            nueY = hitE[1] + detectorPosition[1]
            nueZ = hitE[2] + detectorPosition[2]
            dsNue = math.sqrt((xDcy-nueX)**2 + (yDcy-nueY)**2 + (zDcy-nueZ)**2)
            sNue = sDcy + dsNue
            tNue = tDcy + dsNue*1E9/c
            if ((abs(nueX) < 2.50) and (abs(nueY) < 2.50)):
                eW = eventWeight
                ENueSpectra = math.sqrt(nuePx*nuePx + nuePy*nuePy + nuePz*nuePz)
                hNueDet.Fill(ENueSpectra)
            else:
                eW = 0.0
            nueDetector = particle.particle(runNumber, event, sNue, nueX, nueY, nueZ, nuePx, nuePy, nuePz, tNue, eW, "nue")
            if (self.__history): eH.addParticle("nueDetector", nueDetector)
            if (self._muDcyCount < printLimit): print ("nue at detector")


    def muDcyCount(self):
        return deepcopy(self._muDcyCount)


if __name__ == "__main__" :

    __mainPrint = False

    __history = False

    parser = argparse.ArgumentParser()
    parser.add_argument('--dict', help='Select which run condition dictionary to use. Default: MuRingDcy[.dict]', default='MuRingDcy')
    parser.add_argument('--run', help='Select which run number to use. If no run number is specified consecutive run number is used.', default='0')
    parser.add_argument('--studyname', help='Select which study name to use. If no study name is specified study name from dictionary is used.', default='noName')
    parser.add_argument('--p0', help='Select central pion momentum to be generated. If no central pion momentum is specified pion momentum from dictionary is used.',default='0.')
    parser.add_argument('--Mup0', help='Select central muon momentum to be stored in the ring. If no central muon momentum is specified muon momentum from dictionary is used.',default='0.')
    parser.add_argument('--inputFile', help='Specify root file with input histograms. Assumed to be in 31-Target/v1/ Default: target40.root',default='target40.root')
    args = parser.parse_args()

    StudyDir = os.getenv('StudyDir')
    if args.studyname != "noName":
        StudyName = args.studyname
    else:
        StudyName = os.getenv('StudyName')
    print (" StudyDir is ", StudyDir)
    print (" StudyName is ", StudyName)
    
    controlFile = os.path.join(StudyDir, StudyName, args.dict+".dict")
    print (" controlFile is ", controlFile)
 
    ctrlInst = control.control(controlFile)

    #Initialise central pion and muon momenta
    if float(args.p0) != 0.:
        pionMom = float(args.p0)
    else:
        pionMom = ctrlInst.PPi()

    if float(args.Mup0) != 0.:
        muonMom = float(args.Mup0)
    else:
        muonMom = ctrlInst.PMu()

    nuSTRMCnst = nuSTORMConst.nuSTORMConst()

    tlCmplxLength = nuSTRMCnst.TrfLineCmplxLen()
    tlCmplxAngle = nuSTRMCnst.TrfLineCmplxAng()

    normInst = normalisation( __history, muonMom)
# run number needed because it labels various files
    if int(args.run) != 0:
        ctrlInst.setRunNumber(int(args.run))
    runNumber = ctrlInst.runNumber(True)


#       logfile initialisation
    logging.basicConfig(filename=ctrlInst.logFile(), level=logging.INFO) #, encoding='utf-8'

#       Histogram initialisation
    hm = histoManager.histoManager()
    hBins = 120
    hLower = 0.0
    hUpper = 7.2
    hNumuDet = hm.book("Enumu Energy", hBins, hLower, hUpper)
    hNueDet = hm.book("Enue Energy", hBins, hLower, hUpper)
    hTitle = "EnumuNorm"
    eNumuDataNorm = hm.book(hTitle, hBins, hLower, hUpper)
    hTitle = "EnueNorm"
    eNueDataNorm = hm.book(hTitle, hBins, hLower, hUpper)
#       histograms for rate study
    rBins = 110
    rLower = 1600.0
    rUpper = 3800.0
    numuDetTime = hm.book("numu Detector Time", rBins, rLower, rUpper)
    nueDetTime = hm.book("nue Detector Time", rBins, rLower, rUpper)
    rTitle = "numuRate"
    numuRate = hm.book(rTitle, rBins, rLower, rUpper)
    rTitle = "nueRate"
    nueDataRate = hm.book(rTitle, rBins, rLower, rUpper)

#    hmuDcyPnt = hm.book("mu decay z (m)", 100, 0.0, 300.0)
#    hLifetime = hm.book("life time", 100, 0.0, 1000.0)
#  start message
    print ("========  Normalisation run: start  ======== Version ", normalisation.__version__)
    print()
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    logging.info("========  Normalisation run: start  ======== Version %s, ... %s", normalisation.__version__, dt_string)

#   Set up the print flags
#    __mainPrint = ctrlInst.mainPrnt()
    print("__mainPrint is ", __mainPrint)
    if (__mainPrint):
        print("main print set to true")
    else:
        print("main print set to false")


# set up the processing flags
    tlFlag = ctrlInst.tlFlag()
    psFlag = ctrlInst.psFlag()
    lstFlag = ctrlInst.lstFlag()
    muDcyFlag = ctrlInst.muDcyFlag()
    FlshAtDetFlg = ctrlInst.flashAtDetector()
    PSMuonsFlag = ctrlInst.PSMuons()
    ringMuonsFlag = ctrlInst.ringMuons()
    tEqualsZeroFlag = ctrlInst.tEqualsZero()
    pencilBeamFlag = ctrlInst.pencilBeam()
    pDistInputFlag = ctrlInst.pDistInput()
    psDistInputFlag = ctrlInst.psDistInput()

    print (f"Processing flags -- tlflag: {tlFlag} / psFlag: {psFlag} / lstFlag: {lstFlag} / muDcyFlag: {muDcyFlag} / FlshAtDetFlg: \
        {FlshAtDetFlg} / PSMuonsFlag: {PSMuonsFlag} / ringMuonsFlag: {ringMuonsFlag} / pencilBeamFlag: {pencilBeamFlag} / tEqualsZeroFlag: \
        {tEqualsZeroFlag} / pDistInputFlag: {pDistInputFlag} / psDistInputFlag: {psDistInputFlag}")
    logging.info("Processing flags -- tlflag: %s,  psFlag: %s,  lstFlag: %s,  muDcyFlag: %s, FlshAtDetFlg: %s, PSMuonsFlag: %s, ringMuonsFlag: %s", \
        tlFlag, psFlag, lstFlag, muDcyFlag, FlshAtDetFlg, PSMuonsFlag, ringMuonsFlag)
    logging.info("     tEqualsZero: %s", tEqualsZeroFlag)
    logging.info("     pencilBeam: %s, pDistInput: %s, psDistInput: %s", pencilBeamFlag, pDistInputFlag, psDistInputFlag)

# get constants
    piCnst  = PC.PionConst()
    muCnst  = MC.MuonConst()
    c = piCnst.SoL()
    piMass = piCnst.mass()/1000.0

# initialise run number, number of events to generate, central pion momentum, and event weight
    printLimit = ctrlInst.printLimit()
    crossSection = 50
    nEvents = ctrlInst.nEvents()
    eventWeight = crossSection
    logging.info("Run Number: %s,  nEvents: %s,  pion central momentum: %s,  muon central momentum: %s", runNumber, nEvents, pionMom, muonMom)
    print("Run Number: ", runNumber)
    print("nEvents: ", nEvents)
    print("Pion central momentum: ", pionMom)
    print("Muon central momentum: ", muonMom)

# Get the nuSIM path name and use it to set names for the inputfile and the outputFile
    nuSIMPATH = os.getenv('nuSIMPATH')
    filename  = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
    print ("filename: " , filename)
    rootInputFilename = os.path.join(nuSIMPATH, '31-Target/v1/', args.inputFile)
    rootFilename = os.path.join(StudyDir, StudyName, 'normalisation' + str(ctrlInst.runNumber())+'.root')
    trfCmplxFile = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-TrfLineCmplx-Params-v1.0.csv')
    print ("trfCmplxFile: " , trfCmplxFile)
    print("dictionary: ", controlFile)
    print ("numSIMPATH, filename, rootinputfilename, rootfilename, trfCmplxFile \n", nuSIMPATH, "\n", filename, "\n", rootInputFilename, "\n", rootFilename,
         "\n", trfCmplxFile)
    logging.info("Control File: %s,  ", controlFile  )
    outFilename = rootFilename
    logging.info("Parameters: %s,  \n     transfer line parameters: %s,  \n     histogram input file: %s \n     output file: %s", filename,  trfCmplxFile, rootInputFilename, rootFilename)


# Get machine and run parameters
    geV = int(pionMom)
#   fractional part
    fracE = int((pionMom - geV)*10)
     
    histName = 'histP'+str(geV)+str(fracE)+'GeV'
    histName2Dx = 'histXPS'+str(geV)+str(fracE)+'GeV'
    histName2Dy = 'histYPS'+str(geV)+str(fracE)+'GeV'
    print ("histName is ", histName, ";    histName2Dx is ", histName2Dx, ";    histName2Dy is ", histName2Dy)
    RndmGen = Rndm.RandomGenerator(rootInputFilename,histName,histName2Dx,histName2Dy)
    psLength = nuSTRMCnst.ProdStrghtLen()
    detectorPosZ = nuSTRMCnst.HallWallDist()
    nuTrLnCmplx = nuTrfLineCmplx.nuSTORMTrfLineCmplx(trfCmplxFile)
    print ("Transfer Line Complex object is ", nuTrLnCmplx)
# set up the detector front face
    zPlPos = 230.0
    xPlPos = 0.0
    yPlPos = 0.0
    detectorPosition = [xPlPos, yPlPos, zPlPos]
    fluxPlane = plane.plane(detectorPosition)
#    fluxPlane = plane.plane(psLength, detectorPosZ)
# set up the event history - instantiate
    if (__history): eH = eventHistory.eventHistory()
    if (__history): eH.outFile(outFilename)
    if (__history): eH.cd()
# create the root structure to write to
    if (__history): eH.rootStructure()
# initialise the python arrays to something sensible
    if (__history): eH.makeHistory()

# event loop
for event in range(nEvents):
# generate a pion
    if (event < printLimit): print ("\n ========== Event Start ==========  ")

    pi = piEvtInst.PionEventInstance(pionMom)
# set its values
    tsc = pi.getLclTraceSpaceCoord()
    if __mainPrint:
        print (f"----> main:(1) tsc is {tsc}")
    s = tsc[0]
    zl= tsc[3]

    if (pDistInputFlag):
        pPion = RndmGen.getRandom()
    else:
        pPion = pi.getppiGen()

    if (psDistInputFlag):
        xl, xpl = RndmGen.getRandom2Dx()
        yl, ypl = RndmGen.getRandom2Dy()
    elif (pencilBeamFlag):
        xl = 0.0
        yl = 0.0
        xpl = 0.0
        ypl = 0.0
    else:
        xl = tsc[1]
        yl = tsc[2]
        xpl = tsc[4]
        ypl = tsc[5]

    pzl = np.sqrt(pPion**2/(1+xpl**2+ypl**2))
    pxl = pzl*xpl
    pyl = pzl*ypl

    if (tEqualsZeroFlag):
        t = 0.0
    else:
        t = nuTrLnCmplx.GenerateTime()*1E9
# transform to the global system
    xg, yg, zg, pxg, pyg, pzg = normInst.tltoGlbl(xl, yl, zl, pxl, pyl, pzl)
#  pion at target is a point source but with a momentum spread
    pionTarget = particle.particle(runNumber, event, s, xg, yg, zg, pxg, pyg, pzg, t, eventWeight, "pi+")
    if (__history): eH.addParticle('target', pionTarget)

# generate a pion
    if ((pDistInputFlag) or (psDistInputFlag) or (pencilBeamFlag)):
        pi = piEvtInst.PionEventInstance(particleTar=pionTarget)
        tsc = pi.getLclTraceSpaceCoord()
        if __mainPrint:
            print (f"----> main:(2) tsc is {tsc}")

# get the decay length - in the transfer line, in the production straight - or lost beyond the straight
    lifetime = pi.getLifetime()
    pathLength = lifetime*pPion*c/piMass

# get the muon momentum
    mu = pi.getmu4mmtm()
    eMu  = mu[0]
    pxMu = mu[1][0]
    pyMu = mu[1][1]
    pzMu = mu[1][2]
    pMu  = math.sqrt(pxMu*pxMu + pyMu*pyMu + pzMu*pzMu)
# and decay the muon
    if __mainPrint:
        print (f"----> main: tsc for neutrinoEventInstance is {pi.getTraceSpaceCoord()}")

    nuEvt = nuEvtInst.NeutrinoEventInstance(pMu,pi.getTraceSpaceCoord())

    Absorbed = nuEvt.Absorption(pi.getTraceSpaceCoord(), pi.getmu4mmtm(), pi.getcostheta(), muonMom)

# decay in the transferline
    if ((tlFlag) and (pathLength < tlCmplxLength)):
      normInst.tlDecay()
      if (muDcyFlag):normInst.decayMuons()
    else:
# pion reaches end of transfer line just write out a new pion with altered s and z - all other pions must do this
# the magnets bend the beam so that the local co-ordinates are now the global ones
      se = tlCmplxLength
      ze = 0.0
# t = d/(beta*c)
      Epion = math.sqrt(pPion**2 + piMass**2)
      te = t + 1E9*se*Epion/(c*pPion)
# x local is the same as before.
      pionPS = particle.particle(runNumber, event, se, xl, yl, ze, pxl, pyl, pzl, te, eventWeight, "pi+")
      if (__history): eH.addParticle("productionStraight", pionPS)
# decay beyond the end of the production straight
    if ((lstFlag) and (pathLength > tlCmplxLength + psLength)):
        normInst.beyondPS()

# decay in the production straight
    if ((psFlag) and (pathLength >= tlCmplxLength) and (pathLength <= tlCmplxLength + psLength)):
        normInst.decayPiInPS()
# decay the muons
        if (muDcyFlag): normInst.decayMuons()



#  write to the root structure
    if (__history): eH.fill()
# tell the user what is happening
    if (event < 10):
        print ("event number is ", event)
        print (pi)
        print (tsc)
        print (pionTarget)
    elif ((event <100) and (event%10 ==0)):
        print ("event number is ", event)
    elif ((event <1000) and (event%100 ==0)):
        print ("event number is ", event)
    elif ((event < 10000) and (event%1000 ==0)):
        print ("event number is ", event)
    elif ((event < 100000) and (event%10000 ==0)):
        print ("event number is ", event)
    else:
        if (event%100000 ==0): print ("event number is ", event)

    if (event == nEvents-1):
        print()
        print(normInst.muDcyCount()," neutrinos have been created.")

#   Normalisation code
#   Constants for normalisation
normValidate = True

nPions = nEvents
# 1 GeV: 0.1163    2 Gev: 0.1338    3 GeV: 0.1240    4 GeV: 0.1120    5 GeV: 0.0985    6 GeV: 0.0865    7 GeV: 0.0760    8 GeV: 0.0671
ppp = [0.1163, 0.1338, 0.1240, 0.1120, 0.0985, 0.0865, 0.0760, 0.0671]
#   Get the pion momentum from the input file
ePi = int(pionMom)
pionPerProton = ppp[ePi-1]
if (normValidate):
    print ("Number of events  is  ", nPions)    
    print ("pion momentum is ", ePi, "     pions per proton is ", pionPerProton)
#   Start with the pion per proton error being 5%
pionPerProtonError = 0.00
#   put the pion to neutrino error to 2%
neutrinoPerPionError = 0.00
#   total uncertainty
uncertainty = math.sqrt(pionPerProtonError*pionPerProtonError + neutrinoPerPionError*neutrinoPerPionError)
if (normValidate):
    print ("uncertainty is  ", uncertainty)
#   Muon neutrinos
nNuMu = hNumuDet.GetEntries()
if (normValidate):
    print ("nNuMu ", nNuMu)
nuMuPerPion = nNuMu/nPions
nuMuPerProton = nuMuPerPion*pionPerProton
effectiveProtons = nPions/pionPerProton
runNorm = 1000000/effectiveProtons
nuMuPerMillProton = nuMuPerProton*runNorm
if (normValidate):
    print ("effectiveProtons ", effectiveProtons)
    print ("runNorm ", runNorm)
    print ("nuMuPerProton ", nuMuPerProton)
    print ("nuMuPerMillProton ", nuMuPerMillProton)
xAxMu = hNumuDet.GetXaxis()
xAxMu.SetTitle("Neutrino Energy (GeV)")
yAxMu = hNumuDet.GetYaxis()
yAxMu.SetTitle("Neutrinos")
hNumuDet.SetTitle("#nu_{#mu} Energy at the detector front face")
hNumuDet.SetStats(0)

nuMu = []
nuMuStatEr = []
nuMuTotalPerCent = []
nuMuTotalErr = []
for pnt in range(0,hNumuDet.GetNbinsX()):
  nuMu.append(hNumuDet.GetBinContent(pnt))
  if nuMu[pnt] > 0.0:
    nuMuStatEr.append(math.sqrt(nuMu[pnt])/nuMu[pnt])
  else:
    nuMuStatEr.append(0.0)

if (normValidate):
    print ("nuMu is \n ", nuMu)

for pnt in range(0,hNumuDet.GetNbinsX()):
  eNumuDataNorm.SetBinContent(pnt, nuMu[pnt]*nuMuPerMillProton)
  tE = math.sqrt(nuMuStatEr[pnt]*nuMuStatEr[pnt] + uncertainty*uncertainty)
  nuMuTotalPerCent.append(tE)
  eNumuDataNorm.SetBinError(pnt, nuMuTotalPerCent[pnt]*nuMu[pnt]*nuMuPerMillProton)

if (normValidate):
    print ("eNuMuDataNorm is \n ", eNumuDataNorm)
    eNumuDataNorm.Print("All")  

#   histogram titles
xAxMuNorm = eNumuDataNorm.GetXaxis()
xAxMuNorm.SetTitle("Neutrino Energy (GeV)")
yAxMuNorm = eNumuDataNorm.GetYaxis()
yAxMuNorm.SetTitle("neutrinos per 10^{6} protons")
eNumuDataNorm.SetTitle("#nu_{#mu} Energy at the detector front face")


#   electron neutrinos
nNuE  = hNueDet.GetEntries()
print ("nNuE ", nNuE)
nuEPerPion = nNuE/nPions
nuEPerProton = nuEPerPion*pionPerProton
nuEPerMillProton = nuEPerProton*runNorm
if (normValidate):
    print ("nuEPerPion ", )
    print ("nuEPerProton ", nuEPerProton)
    print ("nuEPerMillProton ", nuEPerMillProton)


print ("nuEPerProton ", nuEPerProton)
xAxE = hNueDet.GetXaxis()
xAxE.SetTitle("Neutrino Energy (GeV)")
yAxE = hNueDet.GetYaxis()
yAxE.SetTitle("Neutrinos")
hNueDet.SetTitle("#nu_{e} Energy at the detector front face")
hNueDet.SetStats(0)

#   Calculate the errors on the histogram
nuE = []
nuEStatEr = []
nuETotalPerCent = []
nuETotalErr = []
for pnt in range(0,hNueDet.GetNbinsX()):
  nuE.append(hNueDet.GetBinContent(pnt))
  if nuE[pnt] > 0.0:
    nuEStatEr.append(math.sqrt(nuE[pnt])/nuE[pnt])
  else:
    nuEStatEr.append(0.0)

for pnt in range(0,hNueDet.GetNbinsX()):
  eNueDataNorm.SetBinContent(pnt, nuE[pnt]*nuEPerMillProton)
  tE = math.sqrt(nuEStatEr[pnt]*nuEStatEr[pnt] + uncertainty*uncertainty)
  nuETotalPerCent.append(tE)
  eNueDataNorm.SetBinError(pnt, nuETotalPerCent[pnt]*nuE[pnt]*nuEPerMillProton)

xAxENorm = eNueDataNorm.GetXaxis()
xAxENorm.SetTitle("Neutrino Energy (GeV)")
yAxENorm = eNueDataNorm.GetYaxis()
yAxENorm.SetTitle("neutrinos per 10^{6} protons")
eNueDataNorm.SetTitle("#nu_{e} Energy at the detector front face")
eNueDataNorm.SetStats(0)
#hNueDet.Scale(nuEPerProton)

#       Calculate the rate.
#   Data Range is (upper limit - lower limit) ns
#   Bin widths is (Data Range)/(number of Bins) in ns
#   Multiply the number of events in a bin by 1/(bin width) to get the rate in Hz
#   This is the number for the number of pions in the run = nPions
#   This implies the number of protons nProtons = nPions
#   Number of Muon neutrinos is 


#   Number of protons in a spill

totalNumu = numuDetTime.GetEntries()
print ("----> Rate <---- totalNumu at detector ", totalNumu)
numProtons = nPions/pionPerProton
print ("----> Rate <---- pions per proton ", pionPerProton, "    protons for the run is ", numProtons)

#   So a burst with that many protons we would get the number of neutrinos seen, but we want a spill
#   of 5*10^12 protons
protonsPerSpill = 5*10**12
print ("----> Rate <---- protons per spill ", protonsPerSpill)
#   Factor to multiple the neutrino count by to give the rate per proton


#   Size of time bin - need to get directly from the plots - this will be in ns
tUpper = numuDetTime.GetXaxis().GetXmax()
tLower = numuDetTime.GetXaxis().GetXmin()
rangeT = rUpper-rLower
tBins = numuDetTime.GetNbinsX()
tBinWidth = rangeT*10**-9/tBins 
print ("----> Rate <---- t bin width ", tBinWidth)

#   Loop over the muon histogram to get values and calculate statistical error
muCount = []
muCountErr = []
muRate = []
muRateErr = []

#   Unpack histogram
for pnt in range(0,tBins):
  muCount.append(numuDetTime.GetBinContent(pnt))
  if muCount[pnt] > 0.0:
    muCountErr.append(math.sqrt(muCount[pnt]))
  else:
    muCountErr.append(0.0)

#   Fill rate histogram - mutliply by protons/spill and divide by protons in the generation
#   and to get the rate divide by the bin width
for pnt in range(0,tBins):
  numuRate.SetBinContent(pnt, muCount[pnt]*protonsPerSpill/(numProtons*tBinWidth))
  numuRate.SetBinError(pnt, muCountErr[pnt]*protonsPerSpill/(numProtons*tBinWidth))

#   Get the total number of neutrinos
nuTotal = numuDetTime.GetEntries()
print ("----> Rate <---- total number of neutrinos ", nuTotal)
#   In 2ms - or 2000ns. Really should get this from the histogram
totalTRange = 2000*10**-9
print ("----> Rate <---- total range ",totalTRange)
rateTotal = nuTotal*protonsPerSpill/(numProtons*totalTRange)
print ("----> Rate <---- total Rate ",rateTotal)
rateErr = math.sqrt(nuTotal)*protonsPerSpill/(numProtons*totalTRange)
perCentErr = rateErr*100.0/rateTotal
print ("----> Rate <----  percentage rate Error ",perCentErr)

#   From flux to event rate.
Avagadro = 6.0*10**23
dVol = 125              # cubic metres
argonDensity = 1400     # 1.4 gm/l is 1400 gm/cu m
moleWt = 40
nNucleii = dVol*argonDensity*Avagadro/moleWt
print ("----> Rate <----  number of nucleii ", nNucleii)
crossSection = 1.0*10**-38
print ("----> Rate <----  muon neutrino - Argon cross-section ", crossSection)

#   Interaction rate
scatterRate = rateTotal*nNucleii*crossSection
print ("----> Rate <----  muon neutrino - scattering Rate ", scatterRate)

#   Interactions per spill
interPerSpill = scatterRate*totalTRange
print ("----> Rate <----  muon neutrino - interactions per Spill ", interPerSpill)




# Write to the root output file and close
if (__history): eH.write()
if (__history): eH.outFileClose()
# Write out histograms
fileName = os.path.join(StudyDir, StudyName + "/Normalplots" + str(ctrlInst.runNumber()) + ".root")
print(f"filename is {fileName}")
print (fileName)
hm.histOutRoot(fileName)
##! Complete:
print()
print("========  Normalisation run : complete  ========")
