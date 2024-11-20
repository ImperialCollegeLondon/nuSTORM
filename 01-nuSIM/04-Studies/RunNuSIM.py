#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model for calculating normalised numbers
========================================

    Assumes that nuSim code is in python path.

    @file RunNuSIM.py - copiedd from RunCmpltSimulation.py version 1.4

    @brief   calculates the normalisation of the data

    @author  Paul Kyberd

    @version     1.0
    @date        05 October 2021
    @author     Paul Kyberd

    Add python logging
    @version     1.1
    @date        07 January 2022
    @author     Paul Kyberd

    Add the event history at the end of the production straight
    @version     1.2
    @date        07 January 2022
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

    change the initialisation of plane to allow the definition of the position in 
    x,y and z
    @version    1.4
    @date       27 July 2020
    @author     Paul Kyberd

    Introduce the genPion class to clean up the implementation of the different
    running conditions
    @version    1.5
    @date       17 October 2024
    @author     Paul Kyberd

    Modification so a pion which decays in the transfer line is not shown at the start of
    the production straight - and nor is anything else
    @version    1.6
    @date       23 October 2024
    @author     Paul Kyberd

    Modification so a pion which decays after the end of the production straight is shown at
    the end of the production straight
    @version    1.7
    @date       23 October 2024
    @author     Paul Kyberd


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
import genPion as genPion

class normalisation:

    __version__ = 1.7
    def __init__(self, muonMom=0):
        self._tlDcyCount = 0
        self._byndPSCount = 0
        self._PSDcyCount = 0
        self._muDcyCount = 0
        self._tlAngle = tlCmplxAngle*math.pi/180.0
        self._sth = math.sin(self._tlAngle)
        self._cth = math.cos(self._tlAngle)
        self._mup0 = muonMom

# transform x,y,z and px,py,pz co-ordinates from the transfer line local co-ordinates to the
# global ones

    def tltoGlbl(self, xl, yl, zl, pxl, pyl, pzl):

        xg = xl*self._cth + zl*self._sth
        yg = yl
        zg = zl*self._cth - xl*self._sth

        pxg = pxl*self._cth + pzl*self._sth
        pyg = pyl
        pzg = pzl*self._cth - pxl*self._sth

        return xg, yg, zg, pxg, pyg, pzg

# Deal with decay in the transfer line
    def tlDecay(self):

      if (self._tlDcyCount < printLimit): print("----------------> tlDecay: starting")

      self._tlDcyCount = self._tlDcyCount + 1
# decays at present just set the productionStraight particle with a weight of zero - but updated s and z (small value of pz
# for tracespace calculation)
      noParticle = particle.particle(runNumber, event, tlCmplxLength, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, "none")
      eH.addParticle("productionStraight", noParticle)
# add a pion decay particle - set the time to the decay lifetime
      if (self._tlDcyCount < printLimit): print ("--------> tlDecay: pi in tlDecay ", pi)
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
      if (self._tlDcyCount < printLimit): print ("--------> tlDecay: about to add a particle")
      eH.addParticle("pionDecay", pionTLDecay)
      if (self._tlDcyCount < printLimit): print ("--------> tlDecay:  pion at pionDecay: TL ", pionTLDecay)
# add the pion flash neutrino
      numu = pi.getnumu4mmtm()
      pxnul = numu[1][0]
      pynul = numu[1][1]
      pznul = numu[1][2]
      xd, yd, zd, pxnu, pynu, pznu = self.tltoGlbl(xdl, ydl, zdl, pxnul, pynul, pznul)
      nuFlashTL = particle.particle(runNumber, event, sd, xd, yd, zd, pxnu, pynu, pznu, td, eventWeight, "numu")
      eH.addParticle("piFlashNu", nuFlashTL)
      if (self._tlDcyCount < printLimit): print ("--------> tlDecay: at piFlashNu: TL")
# add the muon from the pion flash ... not going to track this, weight non zero, but don't track further ?
      mu = pi.getmu4mmtm()
      pxmul = mu[1][0]
      pymul = mu[1][1]
      pzmul = mu[1][2]
      xd, yd, zd, pxmu, pymu, pzmu = self.tltoGlbl(xdl, ydl, zdl, pxmul, pymul, pzmul)
      muonProdTL = particle.particle(runNumber, event, sd, xd, yd, zd, pxmu, pymu, pzmu, td, eventWeight, "mu+")
      eH.addParticle("muonProduction",muonProdTL)
      if (self._tlDcyCount < printLimit): print ("--------> tlDecay: at muonProduction: TL")
# extrapolate the neutrino the the detector plane
      hitMu = fluxPlane.findHitPositionPiFlash(nuFlashTL)
      if (self._tlDcyCount < printLimit): print ("--------> tlDecay: hit position of neutrino ", hitMu)
#  fill the event History
      numuX = hitMu[0]
      numuY = hitMu[1]
      numuZ = hitMu[2]
      dsNumu = math.sqrt((xd-numuX)**2 + (yd-numuY)**2 + (zd-numuZ)**2)
      sNumu = sd + dsNumu
      tNumu = td + dsNumu*1E9/c + t
      if (self._tlDcyCount < printLimit): print ( "--------> tlDecay: sNumu is ", sNumu, "    dsNumu is ", dsNumu, "     tNumu is ", tNumu, "    c is ", c)
      if ((abs(numuX) < 100) and (abs(numuY) < 100)):
          eW = eventWeight
      else:
          eW = 0.0
      numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, pxnu, pynu, pznu, tNumu, eW, "numu")
      eH.addParticle("numuDetector", numuDetector)
      if (self._tlDcyCount < printLimit): print ("--------> tlDecay: numu at detector")

      return
#
# Deal with decay beyond the production straight -------------------------------------------------
#
    def beyondPS(self):

      if (self._byndPSCount < printLimit): print("----------------> beyondPS: starting")

      self._byndPSCount = self._byndPSCount + 1
      # extraoplate the muon to the end of the production straight??
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
      eH.addParticle("prodStraightEnd", piEnd)
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
      if (self._byndPSCount < printLimit): print ("--------> beyondPS: pi in beyondPS: decayLength ", sd)
      pionLostDecay = particle.particle(runNumber, event, sd, xd, yd, zd, pxd, pyd, pzd, td, eventWeight, "none")
      eH.addParticle("pionDecay", pionLostDecay)
      if (self._byndPSCount < printLimit):  print ("--------> beyondPS: at pionDecay lost")
# add the pion flash neutrino ... set everything to zero - including eventWeight
      nuFlashLost = particle.particle(runNumber, event, sd, xd, yd, zd, pxd, pyd, pzd, td, 0.0, "none")
      if (self._byndPSCount < printLimit):  print ("--------> beyondPS: at piFlashNu lost")
      eH.addParticle("piFlashNu", nuFlashLost)
# add the muon from the pion flash ... set everything to zero - including eventWeight
      muonProdLost = particle.particle(runNumber, event, sd, xd, yd, zd, pxd, pyd, pzd, td, 0.0, "none")
      if (self._byndPSCount < printLimit):  print ("--------> beyondPS: at muonProduction lost")
      eH.addParticle("muonProduction",muonProdLost)
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

      if (self._PSDcyCount < printLimit): print("----------------> decayPiInPS: starting")

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
        print ("--------> decayPiInPS: error in the time")
#      hTotal.Fill(td)
#      hLifetime.Fill(piLifetime)
#      hStarttime.Fill(t)
#      hS.Fill(sd)

      eH.addParticle("pionDecay", pionPSDecay)
      if (self._PSDcyCount < printLimit): print ("--------> decayPiInPS: pionDecay in PS")
# add the muon
      mu = pi.getmu4mmtm()
      pxMu = mu[1][0]
      pyMu = mu[1][1]
      pzMu = mu[1][2]
      eMu = mu[0]
      muProd = particle.particle(runNumber, event, sd, xd, yd, zd, pxMu, pyMu, pzMu, td, eventWeight, "mu+")
      eH.addParticle("muonProduction", muProd)
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
          eH.addParticle("prodStraightEnd", muEnd)
      else:
          noParticle = particle.particle(runNumber, event, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0,   "none")
          eH.addParticle("prodStraightEnd", noParticle)

      if (self._PSDcyCount < printLimit): print ("--------> decayPiInPS: muonProduction in production straight")
# add the pion flash neutrino
      numu = pi.getnumu4mmtm()
      pxnu = numu[1][0]
      pynu = numu[1][1]
      pznu = numu[1][2]
      eNu = numu[0]

      nuFlash = particle.particle(runNumber, event, sd, xd, yd, zd, pxnu, pynu, pznu, td, eventWeight, "numu")
      eH.addParticle("piFlashNu", nuFlash)
      if (self._PSDcyCount < printLimit):  print ("--------> decayPiInPS: piFlashNu PS")
# extrapolate the neutrino the the detector plane
      if FlshAtDetFlg:
          hitMu = fluxPlane.findHitPositionPiFlash(nuFlash)
          if (self._PSDcyCount < printLimit): print ("--------> decayPiInPS: hit position of neutrino ", hitMu)
#  fill the event History
          numuX = hitMu[0]
          numuY = hitMu[1]
          numuZ = hitMu[2]
          dsNumu = math.sqrt((xd-numuX)**2 + (yd-numuY)**2 + (zd-numuZ)**2)
          sNumu = sd + dsNumu
          tNumu = td + dsNumu*1E9/c + t
          if (self._PSDcyCount < printLimit): print ( "--------> decayPiInPS: neutrino co-ords at detector sNumu is ", sNumu, "    dsNumu is ", dsNumu, "     tNumu is ", tNumu, "    c is ", c)
          if ((abs(numuX) < 2.5) and (abs(numuY) < 2.5)):
            eW = eventWeight
          else:
            eW = 0.0
          numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, pxnu, pynu, pznu, tNumu, eW, "numu")
          eH.addParticle("numuDetector", numuDetector)
      if (self._PSDcyCount < printLimit): print ("--------> decayPiInPS: numu at detector")

#
# Muon decays -------------------------------------------------
#
    def decayMuons(self):

# ... We can allow any muon which is created in the production straight to decay without worrying about
# ... acceptance

      piTraceSpaceCoord = pi.getTraceSpaceCoord()
      mucostheta = pi.getcostheta()
      mu4mom = pi.getmu4mmtm()
      if (self._muDcyCount < printLimit):
            print ("piTraceSpaceCoord is ", piTraceSpaceCoord)
            print ("muon 4 momentum is ", mu4mom)
            print ('muoncostheta is ' , mucostheta)

      mup0 = self._mup0
      if self._mup0 == 0.:
          print("decayMuons: ERROR: Muon central momentum is 0.! Interrupting simulation...")
          exit()

      Absorbed = nuEvt.Absorption(piTraceSpaceCoord, mu4mom, mucostheta, mup0)
#        if (Absorbed == False): print ("absorbed is ", Absorbed, "  for event ", event)
#  if the muon doesn't make it in the ring acceptance ... it is absorbed ... so for the
#  muon we create a suitable muon - but the other particles are all put to null values so
#  all values bar run and event are put to zero
      if (Absorbed):
        testParticle = particle.particle(runNumber, event, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0,   "none")
        eH.addParticle("muonDecay", testParticle)
        if (self._muDcyCount < printLimit): print ("--------> decayMuons: absorbed")
        if (PSMuonsFlag):
          muTSC = nuEvt.getTraceSpaceCoord()
          sDcy = muTSC[0]
#TL decay
#          if (sDcy < tlCmplxLength):
# ps decay
          if (sDcy < tlCmplxLength+psLength):
            print (f"decayMuons: muon in production: {sDcy}")
# Muon Decay
            xDcy = muTSC[1]
            yDcy = muTSC[2]
            zDcy = muTSC[3]
            pxDcy = nuEvt.getPb()[0]
            pyDcy = nuEvt.getPb()[1]
            pzDcy = nuEvt.getPb()[2]
            tDcy = (pi.getLifetime()+nuEvt.getLifeTime())*1E9 + t
            muDecay = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, pxDcy, pyDcy, pzDcy, tDcy, eventWeight, "mu+")
            eH.addParticle("muonDecay", muDecay)
            self._muDcyCount = self._muDcyCount + 1
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: muDecay is ", muDecay)

# electron production
            e4mmtm = nuEvt.gete4mmtm()
            e3mmtm = e4mmtm[1]
            eProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, e3mmtm[0], e3mmtm[1], e3mmtm[2], tDcy, eventWeight, 'e+')
            eH.addParticle('eProduction', eProd)
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: eProduction is ", eProd)

# numu production
            numu4mmtm = nuEvt.getnumu4mmtm()
            numu3mmtm = numu4mmtm[1]
            numuPx = numu3mmtm[0]
            numuPy = numu3mmtm[1]
            numuPz = numu3mmtm[2]
            numuProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, numu3mmtm[0], numu3mmtm[1], numu3mmtm[2], tDcy, eventWeight, 'numuBar')
            eH.addParticle('numuProduction', numuProd)
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: numuProduction is ", numuProd)

# nue production
            nue4mmtm = nuEvt.getnue4mmtm()
            nue3mmtm = nue4mmtm[1]
            nuePx = nue3mmtm[0]
            nuePy = nue3mmtm[1]
            nuePz = nue3mmtm[2]
            nueProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, nue3mmtm[0], nue3mmtm[1], nue3mmtm[2], tDcy, eventWeight, 'nue')
            eH.addParticle('nueProduction', nueProd)
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: nueProduction is ", nueProd)
# and finally extrapolate to the neutrino detector
#   hitx is x, y, z, R, phi, px, py, pz, E
            hitE,hitMu=fluxPlane.findHitPositionMuEvt(nuEvt)
#  numu extrapolation to the detector
            numuX = hitMu[0]
            numuY = hitMu[1]
            numuZ = hitMu[2]
            dsNumu = math.sqrt((xDcy-numuX)**2 + (yDcy-numuY)**2 + (zDcy-numuZ)**2)
            sNumu = sDcy + dsNumu
            tNumu = tDcy + sNumu/c + t
            if (self._muDcyCount < printLimit): print ( "--------> decayMuons(a): sNumu is ", sNumu, "    dsNumu is ", dsNumu, "     tNumu is ", tNumu)
            if ((abs(numuX) < 2.50) and (abs(numuY) < 2.50)):
                eW = eventWeight
            else:
                eW = 0.0
            numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, numuPx, numuPy, numuPz, tNumu, eW, "numu")
            eH.addParticle("numuDetector", numuDetector)
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: numu at detector")
#  nue extrapolation to the detector
            nueX = hitE[0]
            nueY = hitE[1]
            nueZ = hitE[2]
            dsNue = math.sqrt((xDcy-nueX)**2 + (yDcy-nueY)**2 + (zDcy-nueZ)**2)
            sNue = sDcy + dsNue
            tNue = tDcy + sNue/c + t
            if ((abs(nueX) < 2.50) and (abs(nueY) < 2.50)):
                eW = eventWeight
            else:
                eW = 0.0
            nueDetector = particle.particle(runNumber, event, sNue, nueX, nueY, nueZ, nuePx, nuePy, nuePz, tNue, eW, "nue")
            eH.addParticle("nueDetector", nueDetector)
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: nue at detector")

# finished dealing with muon decays in the production straight
      else:
        if (ringMuonsFlag):
            self._muDcyCount = self._muDcyCount + 1
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: not absorbed")
            muTSC = nuEvt.getTraceSpaceCoord()
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: muTSC is ", muTSC)
# Muon Decay
            sDcy = muTSC[0]
            xDcy = muTSC[1]
            yDcy = muTSC[2]
            zDcy = muTSC[3]
            pxDcy = nuEvt.getPb()[0]
            pyDcy = nuEvt.getPb()[1]
            pzDcy = nuEvt.getPb()[2]
            tDcy = (pi.getLifetime()+nuEvt.getLifeTime())*1E9 + t
            muDecay = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, pxDcy, pyDcy, pzDcy, tDcy, eventWeight, "mu+")
            eH.addParticle("muonDecay", muDecay)
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: muDecay is ", muDecay)

# electron production
            e4mmtm = nuEvt.gete4mmtm()
            e3mmtm = e4mmtm[1]
            eProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, e3mmtm[0], e3mmtm[1], e3mmtm[2], tDcy, eventWeight, 'e+')
            eH.addParticle('eProduction', eProd)
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: eProduction is ", eProd)

# numu production
            numu4mmtm = nuEvt.getnumu4mmtm()
            numu3mmtm = numu4mmtm[1]
            numuPx = numu3mmtm[0]
            numuPy = numu3mmtm[1]
            numuPz = numu3mmtm[2]
            numuProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, numu3mmtm[0], numu3mmtm[1], numu3mmtm[2], tDcy, eventWeight, 'numuBar')
            eH.addParticle('numuProduction', numuProd)
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: numuProduction is ", numuProd)

# nue production
            nue4mmtm = nuEvt.getnue4mmtm()
            nue3mmtm = nue4mmtm[1]
            nuePx = nue3mmtm[0]
            nuePy = nue3mmtm[1]
            nuePz = nue3mmtm[2]
            nueProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, nue3mmtm[0], nue3mmtm[1], nue3mmtm[2], tDcy, eventWeight, 'nue')
            eH.addParticle('nueProduction', nueProd)
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: nueProduction is ", nueProd)

# and finally extrapolate to the neutrino detector
#   hitx is x, y, z, R, phi, px, py, pz, E
            hitE,hitMu=fluxPlane.findHitPositionMuEvt(nuEvt)
#  numu extrapolation to the detector
            numuX = hitMu[0]
            numuY = hitMu[1]
            numuZ = hitMu[2]
            dsNumu = math.sqrt((xDcy-numuX)**2 + (yDcy-numuY)**2 + (zDcy-numuZ)**2)
            sNumu = sDcy + dsNumu
            tNumu = tDcy + sNumu/c + t
            if (self._muDcyCount < printLimit): print ( "--------> decayMuons(b): sNumu is ", sNumu, "    dsNumu is ", dsNumu, "     tNumu is ", tNumu)
            if (self._muDcyCount < printLimit): print ( "--------> decayMuons(b): numuX is ", numuX, "    numuY is ", numuY)
            if ((abs(numuX) < 2.50) and (abs(numuY) < 2.50)):
                eW = eventWeight
            else:
                eW = 0.0
            numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, numuPx, numuPy, numuPz, tNumu, eW, "numu")
            if (self._muDcyCount < printLimit): print ( "--------> decayMuons: weight is ", eW)
            
            eH.addParticle("numuDetector", numuDetector)
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: numu at detector, weight is ", eW)
#  nue extrapolation to the detector
            nueX = hitE[0]
            nueY = hitE[1]
            nueZ = hitE[2]
            dsNue = math.sqrt((xDcy-nueX)**2 + (yDcy-nueY)**2 + (zDcy-nueZ)**2)
            sNue = sDcy + dsNue
            tNue = tDcy + sNue/c + t
            if ((abs(nueX) < 2.50) and (abs(nueY) < 2.50)):
                eW = eventWeight
            else:
                eW = 0.0
            nueDetector = particle.particle(runNumber, event, sNue, nueX, nueY, nueZ, nuePx, nuePy, nuePz, tNue, eW, "nue")
            eH.addParticle("nueDetector", nueDetector)
            if (self._muDcyCount < printLimit): print ("--------> decayMuons: nue at detector")


    def muDcyCount(self):
        return deepcopy(self._muDcyCount)

class validationHistos:
    '''
    Factorise off the code for plotting histograms during generation
        
    @version    1.0
    @date       22 October 2024
    @author     Paul Kyberd

    '''
    def __init__(self, pionMom):
         
        self._hm = histoManager.histoManager()
        self._hTotal = self._hm.book("total time", 100, 0.0, 11000.0)
        self._hLifetime = self._hm.book("life time", 100, 0.0, 1000.0)
        self._hStarttime = self._hm.book("start time", 100, 0.0, 11000.0)
        self._hS = self._hm.book("s (m)", 100, 0.0, 300.0)
#   new for validation
        hxLimLow = 0.0
        hxLimHigh = 0.0
        hyLimLow = 0.0
        hyLimHigh = 0.0
        hpxLimLow = 0.0
        hpxLimHigh = 0.0
        hpyLimLow = 0.0
        hpyLimHigh = 0.0
        self._hTarS = self._hm.book("s target (m)", 100, -5.0, 5.0)
        self._hTarX = self._hm.book("x (m)", 100, hxLimLow, hxLimHigh)
        self._hTarY = self._hm.book("y (m)", 100, hyLimLow, hyLimHigh)
        self._hTarZ = self._hm.book("z (m)", 100, -52.0, -48.0)
        self._hTarPx = self._hm.book("px (MeV)", 100, hpxLimLow, hpxLimHigh)
        self._hTarPy = self._hm.book("py (MeV)", 100, hpyLimLow, hpyLimHigh)
        self._hTarPz = self._hm.book("pz (MeV)", 100, 0.0, 8.0)
        self._hTarP = self._hm.book("p (MeV)", 100, 0.84*pionMom, 1.16*pionMom)
        self._hTarXp = self._hm.book("Xp", 100, hpxLimLow, hpxLimHigh)
        self._hTarYp = self._hm.book("Yp", 100, hpyLimLow, hpyLimHigh)
        self._hTarXY = self._hm.book2("xy", 100, hxLimLow, hxLimHigh, 100, hyLimLow, hyLimHigh)
        self._hTarXXp = self._hm.book2("x v Xp", 100, hxLimLow, hxLimHigh, 100, hpxLimLow, hpxLimHigh)
        self._hTarYYp = self._hm.book2("y v Yp", 100, hyLimLow, hyLimHigh, 100, hpyLimLow, hpyLimHigh)
        self._hTarXpYp = self._hm.book2("Xp v Yp", 100, hpxLimLow, hpxLimHigh, 100, hpyLimLow, hpyLimHigh)
    
    def Fill(self, s, xl, yl, zl, pxl, pyl, pzl, xpl, ypl):
        p = np.sqrt(pxl*pxl + pyl*pyl +pzl*pzl)
        self._hTarS.Fill(s)
        self._hTarX.Fill(xl)
        self._hTarY.Fill(yl)
        self._hTarZ.Fill(zl)
        self._hTarPx.Fill(pxl)
        self._hTarPy.Fill(pyl)
        self._hTarPz.Fill(pzl)
        self._hTarP.Fill(p)
        self._hTarXp.Fill(xpl)
        self._hTarYp.Fill(ypl)
        self._hTarXY.Fill(xl, yl) 
        self._hTarXXp.Fill(xl, xpl) 
        self._hTarYYp.Fill(yl, ypl)
        self._hTarXpYp.Fill(xpl, ypl)

    def output(self, fileName):
       self._hm.histOutRoot(fileName)


def showEvent(event):

    if (event < 10):
        print ("event number is ", event)
#        print (pi)
#        print (tsc)
#        print (pionTarget)
    elif ((event <100) and (event%10 ==0)):
        print ("event number is ", event)
    elif ((event <1000) and (event%100 ==0)):
        print ("event number is ", event)
    elif ((event < 10000) and (event%1000 ==0)):
        print ("event number is ", event)
    elif ((event < 100000) and (event%10000 ==0)):
        print ("event number is ", event)
    else:
        if (event%100000 == 0): print ("event number is ", event)

def getMomenta(): 

    if float(args.p0) != 0.:
        pionMom = float(args.p0)
    else:
        pionMom = ctrlInst.PPi()

    if float(args.Mup0) != 0.:
        muonMom = float(args.Mup0)
    else:
        muonMom = ctrlInst.PMu()

    return pionMom, muonMom

if __name__ == "__main__" :

    __mainPrint = False


    parser = argparse.ArgumentParser()
    parser.add_argument('--dict', help='Select which run condition dictionary to use. Default: MuRingDcy[.dict]', default='MuRingDcy')
    parser.add_argument('--run', help='Select which run number to use. If no run number is specified consecutive run number is used.', default='0')
    parser.add_argument('--studyname', help='Select which study name to use. If no study name is specified study name from dictionary is used.', default='noName')
    parser.add_argument('--p0', help='Select central pion momentum to be generated. If no central pion momentum is specified pion momentum from dictionary is used.',default='0.')
    parser.add_argument('--Mup0', help='Select central muon momentum to be stored in the ring. If no central muon momentum is specified muon momentum from dictionary is used.',default='0.')
    parser.add_argument('--inputFile', help='Specify root file with input histograms. Default: Scratch/target.root',default='Scratch/target.root')
    args = parser.parse_args()

#   Get the study name and the study directory so that the control file can be accesses
    StudyDir = os.getenv('StudyDir')
    if args.studyname != "noName":
        StudyName = args.studyname
    else:
        StudyName = os.getenv('StudyName')
    print (" main: StudyDir is ", StudyDir)
    print (" main: StudyName is ", StudyName)

#   control file - generate name and open
    controlFile = os.path.join(StudyDir, StudyName, args.dict+".dict")
    print (" main: controlFile is ", controlFile)
    ctrlInst = control.control(controlFile)

# run number needed because it labels various files
    if int(args.run) != 0:
        ctrlInst.setRunNumber(int(args.run))
    runNumber = ctrlInst.runNumber(True)

#       logfile initialisation - using the name in the control file
    logger = logging.getLogger("nuSIM")
    logging.basicConfig(filename=ctrlInst.logFile(), level=logging.INFO) #, encoding='utf-8'
    print ("main: log file is ", ctrlInst.logFile())

    #Initialise central pion and muon momenta
    pionMom, muonMom = getMomenta()
    if __mainPrint:
       print ("main: back from getMomenta: pionMom, muonMom: ", pionMom, "   ", muonMom)

    nuSTRMCnst = nuSTORMConst.nuSTORMConst()

    tlCmplxLength = nuSTRMCnst.TrfLineCmplxLen()
    tlCmplxAngle = nuSTRMCnst.TrfLineCmplxAng()

    normInst = normalisation(muonMom)

#       Histogram initialisation
    vH = validationHistos(pionMom)

#  start message
    print ("========  Normalisation run: start  ======== Version ", normalisation.__version__, "\n")
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    logging.info("========  Normalisation run: start  ======== Version %s, ... %s", normalisation.__version__, dt_string)
#   Set up counters
    piDecayInTL = 0         # Decays in transfer line
    piDecaybeyondPS = 0     # Decays beyond production straight
    piDecayInPS = 0         # Decays in production straight
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
    generator = ctrlInst.genType()
    logging.info("     generator is : %s", generator)

    print (f"Processing flags -- tlflag: {tlFlag} / psFlag: {psFlag} / lstFlag: {lstFlag} / muDcyFlag: {muDcyFlag} / FlshAtDetFlg: \
        {FlshAtDetFlg} / PSMuonsFlag: {PSMuonsFlag} / ringMuonsFlag: {ringMuonsFlag} / pencilBeamFlag: {pencilBeamFlag} / tEqualsZeroFlag: \
        {tEqualsZeroFlag} / pDistInputFlag: {pDistInputFlag} / psDistInputFlag: {psDistInputFlag} / generator: {generator}")
    logging.info("Processing flags -- tlflag: %s,  psFlag: %s,  lstFlag: %s,  muDcyFlag: %s, FlshAtDetFlg: %s, PSMuonsFlag: %s, ringMuonsFlag: %s", \
        tlFlag, psFlag, lstFlag, muDcyFlag, FlshAtDetFlg, PSMuonsFlag, ringMuonsFlag)
    logging.info("     tEqualsZero: %s", tEqualsZeroFlag)
    logging.info("     pencilBeam: %s, pDistInput: %s, psDistInput: %s", pencilBeamFlag, pDistInputFlag, psDistInputFlag)

# get constants
    piCnst  = PC.PionConst()
    c = piCnst.SoL()
    piMass = piCnst.mass()/1000.0

# initialise run number, number of events to generate, central pion momentum, and event weight
    printLimit = ctrlInst.printLimit()
    crossSection = 50
    nEvents = ctrlInst.nEvents()
    logging.info("Run Number: %s,  nEvents: %s,  pion central momentum: %s,  muon central momentum: %s", runNumber, nEvents, pionMom, muonMom)
    print("Run Number: ", runNumber, "\n nEvents: ", nEvents, "\nPion central momentum: ", pionMom, "\nMuon central momentum: ", muonMom)

# Get the nuSIM path name and use it to set names for the inputfile and the outputFile
    nuSIMPATH = os.getenv('nuSIMPATH')
    filename  = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
    mcInputFile = os.path.join(nuSIMPATH, args.inputFile)       # This file should be for the root or KDE input file as appropriate
    rootFilename = os.path.join(StudyDir, StudyName, 'normalisation' + str(ctrlInst.runNumber())+'.root')
    trfCmplxFile = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-TrfLineCmplx-Params-v1.0.csv')
    print ("numSIMPATH, filename, mcInputfile, rootfilename, trfCmplxFile \n", nuSIMPATH, "\n", filename, "\n", mcInputFile, "\n", rootFilename,
         "\n", trfCmplxFile)
    logging.info("Control File: %s,  ", controlFile  )
    outFilename = rootFilename
    logging.info("Parameters: %s,  \n     transfer line parameters: %s,  \n     histogram input file: %s \n     output file: %s", filename,  trfCmplxFile, mcInputFile, rootFilename)

#   instantiate genPion
    gp = genPion.genPion(pionMom, generator, nEvents, mcInputFile)

# Get machine and run parameters
    psLength = nuSTRMCnst.ProdStrghtLen()
    detectorPosZ = nuSTRMCnst.HallWallDist()
    nuTrLnCmplx = nuTrfLineCmplx.nuSTORMTrfLineCmplx(trfCmplxFile)
# set up the detector front face
    zPlPos = 50.0
    xPlPos = 0.0
    yPlPos = 0.0
    position = [xPlPos, yPlPos, zPlPos]
    fluxPlane = plane.plane(position)
#    fluxPlane = plane.plane(psLength, detectorPosZ)
# set up the event history - instantiate
    eH = eventHistory.eventHistory()
    eH.outFile(outFilename)
    eH.cd()
# create the root structure to write to
    eH.rootStructure()
# initialise the python arrays to something sensible
    eH.makeHistory()

# event loop
for event in range(nEvents):
    eventWeight = crossSection
    if __mainPrint:
       print ("\n\n\n\n")
#   Get decay time of the pion
    if (tEqualsZeroFlag):
        t = 0.0
    else:
        t = nuTrLnCmplx.GenerateTime()*1E9

#   genPion - run whatever happens
    pionStart, p = gp.newParticle()
    s = pionStart.s()
    xl = pionStart.x()
    yl = pionStart.y()
    zl = -50.0
    xpl = pionStart.xp()
    ypl = pionStart.yp()
    pzl = np.sqrt(p**2/(1+xpl**2+ypl**2))
    pxl = pzl*xpl
    pyl = pzl*ypl
    pPion = p

#   update validation plots
    vH.Fill(s, xl, yl, zl, pxl, pyl, pzl, xpl, ypl)

# transform to the global system
    xg, yg, zg, pxg, pyg, pzg = normInst.tltoGlbl(xl, yl, zl, pxl, pyl, pzl)

#  Store generated pion in global co-ordinates
    pionTargetGbl = particle.particle(runNumber, event, s, xg, yg, zg, pxg, pyg, pzg, t, eventWeight, "pi+")
    eH.addParticle('target', pionTargetGbl)   
    if __mainPrint:
        print (f"----> main:(2) pionTarget is {pionTargetGbl}")
#   Create a generated pion in local co-ordinates
    pionTarget = particle.particle(runNumber, event, s, xl, yl, zl, pxl, pyl, pzl, t, eventWeight, "pi+")

# generate a pion
    pi = piEvtInst.PionEventInstance(particleTarLocal=pionTarget)
    tsc = pi.getLclTraceSpaceCoord()
    if __mainPrint:
        print (f"----> main:(2) pi is {pi}")
        print (f"----> main:(2) tsc is {tsc}")

# get the decay length - in the transfer line, in the production straight - or lost beyond the straight
    lifetime = pi.getLifetime()
    pathLength = lifetime*pPion*c/piMass
    if __mainPrint:
        print (f"----> main: pathLength is  {pathLength}")

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

#   Now find out where it decays and act appropriately
# decay in the transferline 
    if (pathLength < tlCmplxLength):
        if __mainPrint:
            print (f"----> main: pion decay in transfer lines")
        eventWeight = 0
        piDecayInTL = piDecayInTL + 1
#  only follow if transfer line decay is enabled       
        if (tlFlag):
            normInst.tlDecay()
            if (muDcyFlag):normInst.decayMuons()

#   Decays beyond the end of the transfer line 
    if (pathLength >= tlCmplxLength):
  
# pion reaches end of transfer line just write out a new pion with altered s and z - all other pions must do this
# the magnets bend the beam so that the local co-ordinates are now the global ones
        se = tlCmplxLength
        ze = 0.0
# t = d/(beta*c)
        Epion = math.sqrt(pPion**2 + piMass**2)
        te = t + 1E9*se*Epion/(c*pPion)
# x local is the same as before.
        pionPS = particle.particle(runNumber, event, se, xl, yl, ze, pxl, pyl, pzl, te, eventWeight, "pi+")
        eH.addParticle("productionStraight", pionPS)
        if __mainPrint:
            print (f"----> main: pion at start of production straight {pionPS}")
 
# decay beyond the end of the production straight
    if (pathLength > tlCmplxLength + psLength):
        piDecaybeyondPS = piDecaybeyondPS + 1 
        if __mainPrint:
            print (f"----> main: pion decays beyond end of production straight {pionPS}")
# fill in the end of the production straight even if the rest of the event is not being tracked
        ze1 = 180.0
        se1 = 230.0
        te1 = te + 1E9*ze*Epion/(c*pPion)
        pionPSend = particle.particle(runNumber, event, se1, xl, yl, ze1, pxl, pyl, pzl, te1, eventWeight, "pi+")
        eH.addParticle("prodStraightEnd", pionPSend)
        if (lstFlag):
            normInst.beyondPS()

# decay in the production straight
    if ((psFlag) and (pathLength >= tlCmplxLength) and (pathLength <= tlCmplxLength + psLength)):
        piDecayInPS = piDecayInPS + 1
        normInst.decayPiInPS()
# decay the muons
#   Decays beyond the end of the transfer line 
        if (muDcyFlag): normInst.decayMuons()
        if __mainPrint:
            print (f"----> main: muon decays in the production straight {pionPS}")

#  write to the root structure
    eH.fill()
# tell the user what is happening
    showEvent(event)
## --------- END of EVENT LOOP -----------   


print("\n",normInst.muDcyCount()," neutrinos have been created.")
print("Pions decay in transfer line ", piDecayInTL)
logging.info("Pions decay in transfer line: %i ", piDecayInTL)
print("Pions decay beyond end of production straight ", piDecaybeyondPS)
logging.info("Pions decay beyond end of production straight: %i ", piDecaybeyondPS)
print("Pions decay in production straight ", piDecayInPS)
logging.info("Pions decay in production straight: %i ", piDecayInPS)

# Write to the root output file and close
eH.write()
eH.outFileClose()
# Write out histograms
fileName = os.path.join(StudyDir, StudyName + "/Normalplots" + str(ctrlInst.runNumber()) + ".root")
print(f"filename is {fileName}")
vH.output(fileName)
##! Complete:
print("\n========  Normalisation run : complete  ========")
 