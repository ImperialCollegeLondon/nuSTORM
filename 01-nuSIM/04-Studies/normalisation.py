#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model for calculating normalised numbers
========================================

    Assumes that nuSim code is in python path.

    @file normalisation.py

    @brief   calculates the normalisation of the data

    @author  Paul Kyberd

    Add python logging
    @version     1.1
    @date        07 January 2021


    @version     1.0
    @date        05 October 2021

    @author      Marvin Pfaff
    @version     1.1
    @date        21 January 2022

"""
#
# Class to do the calculation of the event rate normalisation
#
import os, sys
from datetime import datetime
import numpy as np
import math as math
import logging
import PionConst as PC
import nuSTORMConst
import control
import histoManager
import nuSTORMPrdStrght as nuPrdStrt
import nuSTORMTrfLineCmplx as nuTrf
import PionEventInstance as piEvtInst
import NeutrinoEventInstance as nuEvtInst
import RandomGenerator as Rndm
import plane as plane
import particle as particle
import eventHistory as eventHistory

class normalisation:

    __version__ = 1.1
    def __init__(self):
        self._tlDcyCount = 0
        self._byndPSCount = 0
        self._PSDcyCount = 0
        self._muDcyCount = 0
        self._tlAngle = 8.5*math.pi/180.0       # should put the angle in the parameter file
        self._sth = math.sin(self._tlAngle)
        self._cth = math.cos(self._tlAngle)

# transform x,y,z and px,py,pz co-ordinates from the transfer line local co-ordinates to the
# global ones

    def tltoGlbl(self, xl, yl, zl, pxl, pyl, pzl):

        xg = xl + zl*self._sth
        yg = yl
        zg = zl*self._cth

        pxg = pxl*self._cth + pzl*self._sth
        pyg = pyl
        pzg = pzl*self._cth - pxl*self._sth

        return xg, yg, zg, pxg, pyg, pzg

# Deal with decay in the transfer line
    def tlDecay(self):

      self._tlDcyCount = self._tlDcyCount + 1
# decays at present just set the productionStraight particle with a weight of zero - but updated s and z (small value of pz
# for tracespace calculation)
      noParticle = particle.particle(runNumber, event, tlCmplxLength, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, "none")
      eH.addParticle("productionStraight", noParticle)
# add a pion decay particle - set the time to the decay lifetime
      if (self._tlDcyCount < 5): print ("pi in tlDecay ", pi)
      tsc = pi.getTraceSpaceCoord()
      sd = tsc[0]
      xd = tsc[1]
      yd = tsc[2]
      zd = tsc[3]
      xpd = tsc[4]
      ypd = tsc[5]
      pPion = pi.getppiGen()
      pxd = pPion*xpd
      pyd = pPion*ypd
      pzd = np.sqrt(pPion*pPion - pxd**2 - pyd**2)
      pi.getLifetime()
      td = lifetime*1E9 + t
      #zdl = zdl - tlCmplxLength                                                     ##WHAT IS THIS DOING?????
# Now we need to move to the global co-ordinates from the transfer line local
      #xd, yd, zd, pxd, pyd, pzd = self.tltoGlbl(xdl, ydl, zdl, pxdl, pydl, pzdl)
      pionTLDecay = particle.particle(runNumber, event, sd, xd, yd, zd, pxd, pyd, pzd, td, eventWeight, "pi+")
      eH.addParticle("pionDecay", pionTLDecay)
      if (self._tlDcyCount < 5): print (" pion at pionDecay: TL ", pionTLDecay)
# add the pion flash neutrino
      numu = pi.getnumu4mmtm()
      pxnu = numu[1][0]
      pynu = numu[1][1]
      pznu = numu[1][2]
      #xd, yd, zd, pxnu, pynu, pznu = self.tltoGlbl(xdl, ydl, zdl, pxnul, pynul, pznul)
      nuFlashTL = particle.particle(runNumber, event, sd, xd, yd, zd, pxnu, pynu, pznu, td, eventWeight, "numu")
      eH.addParticle("piFlashNu", nuFlashTL)
      if (self._tlDcyCount < 5): print ("at piFlashNu: TL")
# add the muon from the pion flash ... not going to track this, weight non zero, but don't track further ?
      mu = pi.getmu4mmtm()
      pxmu = mu[1][0]
      pymu = mu[1][1]
      pzmu = mu[1][2]
      #xd, yd, zd, pxmu, pymu, pzmu = self.tltoGlbl(xdl, ydl, zdl, pxmul, pymul, pzmul)
      muonProdTL = particle.particle(runNumber, event, sd, xd, yd, zd, pxmu, pymu, pzmu, td, eventWeight, "mu+")
      eH.addParticle("muonProduction",muonProdTL)
      if (self._tlDcyCount < 5): print ("at muonProduction: TL")
# extrapolate the neutrino the the detector plane
      hitMu = fluxPlane.findHitPositionPiFlash(nuFlashTL)
      if (self._tlDcyCount < 5): print ("hit position of neutrino ", hitMu)
#  fill the event History
      numuX = hitMu[0]
      numuY = hitMu[1]
      numuZ = hitMu[2]
      dsNumu = math.sqrt((xd-numuX)**2 + (yd-numuY)**2 + (zd-numuZ)**2)
      sNumu = sd + dsNumu
      tNumu = td + dsNumu*1E9/c + t
      if (self._tlDcyCount < 5): print ( "sNumu is ", sNumu, "    dsNumu is ", dsNumu, "     tNumu is ", tNumu, "    c is ", c)
      if ((abs(numuX) < 40.0) and (abs(numuY) < 40.0)):
          eW = eventWeight
      else:
          eW = 0.0
      numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, pxnu, pynu, pznu, tNumu, eW, "numu")
      eH.addParticle("numuDetector", numuDetector)
      if (self._tlDcyCount < 5): print ("numu at detector")

      return
#
# Deal with decay beyond the production straight -------------------------------------------------
#
    def beyondPS(self):

      self._byndPSCount = self._byndPSCount + 1
#  add a pion decay particle - set the time to the decay lifetime and the s to the pathlength, eventweight to full
# x,y,z to 0.0 and px,py to 0.0, pz to 0.01 so constructor does not
      tsc = pi.getTraceSpaceCoord()
      sd = tsc[0]
      xd = tsc[1]
      yd = tsc[2]
      zd = tsc[3]
      xpd = tsc[4]
      ypd = tsc[5]
      pPion = pi.getppiGen()
      pxd = pPion*xpd
      pyd = pPion*ypd
      pzd = np.sqrt(pPion*pPion - pxdl**2 - pydl**2)
      td = pi.getLifetime()*1E9 + t
      if (self._byndPSCount < 5): print ("pi in beyondPS: decayLength ", sd)
      pionLostDecay = particle.particle(runNumber, event, sd, xd, yd, zd, pxd, pyd, pzd, td, eventWeight, "pi+")
      print("beyond:PS about to add a pion")
      eH.addParticle("pionDecay", pionLostDecay)
      if (self._byndPSCount < 5):  print ("at pionDecay lost")
# add the pion flash neutrino ... set everything to zero - including eventWeight
      pxDaught = 0.0
      pyDaught = 0.0
      pzDaught = 0.01
      nuFlashLost = particle.particle(runNumber, event, sd, xd, yd, zd, pxDaught, pyDaught, pzDaught, td, 0.0, "numu")
      if (self._byndPSCount < 5):  print ("at piFlashNu lost")
      eH.addParticle("piFlashNu", nuFlashLost)
# add the muon from the pion flash ... set everything to zero - including eventWeight
      muonProdLost = particle.particle(runNumber, event, sd, xd, yd, zd, pxDaught, pyDaught, pzDaught, td, 0.0, "mu+")
      if (self._byndPSCount < 5):  print ("at muonProduction lost")
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

      self._PSDcyCount = self._PSDcyCount + 1
      tsc = pi.getTraceSpaceCoord()
      sd = tsc[0]
      xd = tsc[1]
      yd = tsc[2]
      zd = tsc[3]
      xpd = tsc[4]
      ypd = tsc[5]
      pi.getLifetime()
      td = lifetime*1E9 + t
      pPion =pi.getppiGen()
#  need the lifetime in the nuStorm frame in ns
      piLifetime = pi.getLifetime()*1E9*Epion/(piMass)
      td = piLifetime + t
      zd = sd - tlCmplxLength
      pxd = pPion*xpd
      pyd = pPion*ypd
      pzd = np.sqrt(pPion*pPion - pxd**2 - pyd**2)
      pionPSDecay = particle.particle(runNumber, event, sd, xd, yd, zd, pxd, pyd, pzd, td, eventWeight, "pi+")
      if (td < 150.0): print (f"error in the time {td}")
      hTotal.Fill(td)
      hLifetime.Fill(piLifetime)
      hStarttime.Fill(t)
      hS.Fill(sd)
      eH.addParticle("pionDecay", pionPSDecay)
      if (self._PSDcyCount < 5): print ("pionDecay in PS")
# add the muon
      mu = pi.getmu4mmtm()
      pxMu = mu[1][0]
      pyMu = mu[1][1]
      pzMu = mu[1][2]
      eMu = mu[0]
      muProd = particle.particle(runNumber, event, sd, xd, yd, zd, pxMu, pyMu, pzMu, td, eventWeight, "mu+")
      eH.addParticle("muonProduction", muProd)
      if (self._PSDcyCount < 5): print ("muonProduction in production straight")
# add the pion flash neutrino
      numu = pi.getnumu4mmtm()
      pxnu = numu[1][0]
      pynu = numu[1][1]
      pznu = numu[1][2]
      nuFlash = particle.particle(runNumber, event, sd, xd, yd, zd, pxnu, pynu, pznu, td, eventWeight, "numu")
      eH.addParticle("piFlashNu", nuFlash)
      if (self._PSDcyCount < 5):  print ("piFlashNu PS")
# extrapolate the neutrino the the detector plane
      if FlshAtDetFlg:
          hitMu = fluxPlane.findHitPositionPiFlash(nuFlash)
          if (self._PSDcyCount < 5): print ("hit position of neutrino ", hitMu)
#  fill the event History
          numuX = hitMu[0]
          numuY = hitMu[1]
          numuZ = hitMu[2]
          dsNumu = math.sqrt((xd-numuX)**2 + (yd-numuY)**2 + (zd-numuZ)**2)
          sNumu = sd + dsNumu
          tNumu = td + dsNumu*1E9/c + t
          if (self._PSDcyCount < 5): print ( "sNumu is ", sNumu, "    dsNumu is ", dsNumu, "     tNumu is ", tNumu, "    c is ", c)
          if ((abs(numuX) < 10.0) and (abs(numuY) < 10.0)):
              eW = eventWeight
          else:
              eW = 0.0
          numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, pxnu, pynu, pznu, tNumu, eW, "numu")
          eH.addParticle("numuDetector", numuDetector)
      if (self._PSDcyCount < 5): print ("numu at detector")

#
# Muon decays -------------------------------------------------
#
    def decayMuons(self):

# get the muon momentum
      mu = pi.getmu4mmtm()
      pxMu = mu[1][0]
      pyMu = mu[1][1]
      pzMu = mu[1][2]
      pMu = math.sqrt(pxMu*pxMu + pyMu*pyMu + pzMu*pzMu)
# and decay the muon
      nuEvt = nuEvtInst.NeutrinoEventInstance(pMu)

# ... We can allow any muon which is created in the production straight to decay without worrying about
# ... acceptance

      piTraceSpaceCoord = pi.getTraceSpaceCoord()
      mucostheta = pi.getcostheta()
      mu4mom = pi.getmu4mmtm()
      if (self._muDcyCount < 5):
            print ("piTraceSpaceCoord is ", piTraceSpaceCoord)
            print ("muon 4 momentum is ", mu4mom)
            print ('muoncostheta is ' , mucostheta)

      Absorbed = nuEvt.Absorption(piTraceSpaceCoord, mu4mom, mucostheta)
#        if (Absorbed == False): print ("absorbed is ", Absorbed, "  for event ", event)
#  if the muon doesn't make it in the ring acceptance ... it is absorbed ... so for the
#  muon we create a suitable muon - but the other particles are all put to null values so
#  all values bar run and event are put to zero

      if (Absorbed):
        testParticle = particle.particle(runNumber, event, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0,   "mu+")
        eH.addParticle("muonDecay", testParticle)
        if (self._muDcyCount < 5): print ("absorbed")
# here we deal with muons which decay before the end of the production straight and so we do not worry about absorption
        if (PSMuonsFlag):
          muTSC = nuEvt.getTraceSpaceCoord()
          sDcy = muTSC[0]
          if (sDcy < 230.0):
            print (f"muon in production: {sDcy}")
# Muon Decay
            sDcy = muTSC[0]
            xDcy = muTSC[1]
            yDcy = muTSC[2]
            zDcy = muTSC[3]
            pzDcy = nuEvt.getpmu()
            pxDcy = muTSC[4]*pzDcy
            pyDcy = muTSC[5]*pzDcy
            tDcy = (pi.getLifetime()+nuEvt.getLifeTime())*1E9 + t
            muDecay = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, pxDcy, pyDcy, pzDcy, tDcy, eventWeight, "mu+")
            eH.addParticle("muonDecay", muDecay)
            if (self._muDcyCount < 5): print ("muDecay is ", muDecay)

# electron production
            e4mmtm = nuEvt.gete4mmtm()
            e3mmtm = e4mmtm[1]
            eProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, e3mmtm[0], e3mmtm[1], e3mmtm[2], tDcy, eventWeight, 'e+')
            eH.addParticle('eProduction', eProd)
            if (self._muDcyCount < 5): print ("eProduction is ", eProd)
  
# numu production
            numu4mmtm = nuEvt.getnumu4mmtm()
            numu3mmtm = numu4mmtm[1]
            numuPx = numu3mmtm[0]
            numuPy = numu3mmtm[1]
            numuPz = numu3mmtm[2]
            numuProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, numu3mmtm[0], numu3mmtm[1], numu3mmtm[2], tDcy, eventWeight, 'numuBar')
            eH.addParticle('numuProduction', numuProd)
            if (self._muDcyCount < 5): print ("numuProduction is ", numuProd)
 
# nue production
            nue4mmtm = nuEvt.getnue4mmtm()
            nue3mmtm = nue4mmtm[1]
            nuePx = nue3mmtm[0]
            nuePy = nue3mmtm[1]
            nuePz = nue3mmtm[2]
            nueProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, nue3mmtm[0], nue3mmtm[1], nue3mmtm[2], tDcy, eventWeight, 'nue')
            eH.addParticle('nueProduction', nueProd)
            if (self._muDcyCount < 5): print ("nueProduction is ", nueProd)

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
            if (self._muDcyCount < 5): print ( "sNumu is ", sNumu, "    dsNumu is ", dsNumu, "     tNumu is ", tNumu)
            if ((abs(numuX) < 2.50) and (abs(numuY) < 2.50)):
                eW = eventWeight
            else:
                eW = 0.0
            numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, numuPx, numuPy, numuPz, tNumu, eW, "numu")
            eH.addParticle("numuDetector", numuDetector)
            if (self._muDcyCount < 5): print ("numu at detector")
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
            if (self._muDcyCount < 5): print ("nue at detector")

# finished dealing with muon decays in the production straight
      else:
        if (ringMuonsFlag):
            print (f" running ring muons")
            self._muDcyCount = self._muDcyCount + 1
            if (self._muDcyCount < 5): print ("not absorbed")
            muTSC = nuEvt.getTraceSpaceCoord()
            if (self._muDcyCount < 5): print ("muTSC is ", muTSC)
# Muon Decay
            sDcy = muTSC[0]
            xDcy = muTSC[1]
            yDcy = muTSC[2]
            zDcy = muTSC[3]
            pzDcy = nuEvt.getpmu()
            pxDcy = muTSC[4]*pzDcy
            pyDcy = muTSC[5]*pzDcy
            tDcy = (pi.getLifetime()+nuEvt.getLifeTime())*1E9 + t
            muDecay = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, pxDcy, pyDcy, pzDcy, tDcy, eventWeight, "mu+")
            eH.addParticle("muonDecay", muDecay)
            if (self._muDcyCount < 5): print ("muDecay is ", muDecay)

# electron production
            e4mmtm = nuEvt.gete4mmtm()
            e3mmtm = e4mmtm[1]
            eProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, e3mmtm[0], e3mmtm[1], e3mmtm[2], tDcy, eventWeight, 'e+')
            eH.addParticle('eProduction', eProd)
            if (self._muDcyCount < 5): print ("eProduction is ", eProd)

# numu production
            numu4mmtm = nuEvt.getnumu4mmtm()
            numu3mmtm = numu4mmtm[1]
            numuPx = numu3mmtm[0]
            numuPy = numu3mmtm[1]
            numuPz = numu3mmtm[2]
            numuProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, numu3mmtm[0], numu3mmtm[1], numu3mmtm[2], tDcy, eventWeight, 'numuBar')
            eH.addParticle('numuProduction', numuProd)
            if (self._muDcyCount < 5): print ("numuProduction is ", numuProd)

# nue production
            nue4mmtm = nuEvt.getnue4mmtm()
            nue3mmtm = nue4mmtm[1]
            nuePx = nue3mmtm[0]
            nuePy = nue3mmtm[1]
            nuePz = nue3mmtm[2]
            nueProd = particle.particle(runNumber, event, sDcy, xDcy, yDcy, zDcy, nue3mmtm[0], nue3mmtm[1], nue3mmtm[2], tDcy, eventWeight, 'nue')
            eH.addParticle('nueProduction', nueProd)
            if (self._muDcyCount < 5): print ("nueProduction is ", nueProd)

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
            if (self._muDcyCount < 5): print ( "sNumu is ", sNumu, "    dsNumu is ", dsNumu, "     tNumu is ", tNumu)
            if ((abs(numuX) < 10.0) and (abs(numuY) < 10.0)):
                eW = eventWeight
            else:
                eW = 0.0
            numuDetector = particle.particle(runNumber, event, sNumu, numuX, numuY, numuZ, numuPx, numuPy, numuPz, tNumu, eW, "numu")
            eH.addParticle("numuDetector", numuDetector)
            if (self._muDcyCount < 5): print ("numu at detector")
#  nue extrapolation to the detector
            nueX = hitE[0]
            nueY = hitE[1]
            nueZ = hitE[2]
            dsNue = math.sqrt((xDcy-nueX)**2 + (yDcy-nueY)**2 + (zDcy-nueZ)**2)
            sNue = sDcy + dsNue
            tNue = tDcy + sNue/c + t
            if ((abs(nueX) < 10.0) and (abs(nueY) < 10.0)):
                eW = eventWeight
            else:
                eW = 0.0
            nueDetector = particle.particle(runNumber, event, sNue, nueX, nueY, nueZ, nuePx, nuePy, nuePz, tNue, eW, "nue")
            eH.addParticle("nueDetector", nueDetector)
            
            if (self._muDcyCount < 5): print ("nue at detector")


if __name__ == "__main__" :


    controlFile = "101-Studies/pencilValidation/PSMuDcyControlFile.dict"
#    controlFile = "101-Studies/pencilValidation/PSMuRingDcy.dict"
#    controlFile = "101-Studies/pencilValidation/PSnoMuDcyControlFile.dict"
#    controlFile = "101-Studies/pencilValidation/PSPiFlash.dict"
#    controlFile = "101-Studies/pencilValidation/TLPiFlash.dict"

    ctrlInst = control.control(controlFile)
    normInst = normalisation()
# run number needed because it labels various files
    runNumber = ctrlInst.runNumber(True)

#       logfile initialisation
    logging.basicConfig(filename=ctrlInst.logFile(), encoding='utf-8', level=logging.INFO)

#       Histogram initialisation
    hm = histoManager.histoManager()
    hTotal = hm.book("total time", 100, 0.0, 11000.0)
    hLifetime = hm.book("life time", 100, 0.0, 1000.0)
    hStarttime = hm.book("start time", 100, 0.0, 11000.0)
    hS = hm.book("s (m)", 100, 0.0, 300.0)
#  start message
    print ("========  Normalisation run: start  ======== Version ", normalisation.__version__)
    print()
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    logging.info("========  Normalisation run: start  ======== Version %s, ... %s", normalisation.__version__, dt_string)
    logging.info("Control file %s", controlFile)
    logging.info("\n   Description ... %s ",  ctrlInst.description())


# set up the processing flags
    tlFlag = ctrlInst.tlFlag()
    psFlag = ctrlInst.psFlag()
    lstFlag = ctrlInst.lstFlag()
    muDcyFlag = ctrlInst.muDcyFlag()
    FlshAtDetFlg = ctrlInst.flashAtDetector()
    PSMuonsFlag = ctrlInst.PSMuons()
    ringMuonsFlag = ctrlInst.ringMuons()
    print (f"Processing flags -- tlflag: {tlFlag} / psFlag: {psFlag} / lstFlag: {lstFlag} / muDcyFlag: {muDcyFlag} / FlshAtDetFlg: \
        {FlshAtDetFlg} / PSMuonsFlag: {PSMuonsFlag} / ringMuonsFlag {ringMuonsFlag}")
    logging.info("Processing flags -- tlflag: %s,  psFlag: %s,  lstFlag: %s,  muDcyFlag: %s, FlshAtDetFlg: %s, PSMuonsFlag: %s, ringMuonsFlag: %s", \
        tlFlag, psFlag, lstFlag, muDcyFlag, FlshAtDetFlg, PSMuonsFlag, ringMuonsFlag)

# get constants
    piCnst  = PC.PionConst()
    c = piCnst.SoL()
    piMass = piCnst.mass()/1000.0

# initialise run number, number of events to generate, central pion momentum, and event weight
    runNumber = 105
    pionMom = 5.0
    crossSection = 50
    nEvents = ctrlInst.nEvents()
    eventWeight = crossSection
    logging.info("Run Number: %s,  nEvents: %s,  pion central momentum: %s", runNumber, nEvents, pionMom)

# Get the nuSIM path name and use it to set names for the inputfile and the outputFile
    nuSIMPATH = os.getenv('nuSIMPATH')
    filename  = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
    rootInputFilename = os.path.join(nuSIMPATH, 'Scratch/plots_endProdStrght.root')
    rootFilename = os.path.join(nuSIMPATH, 'Scratch/normalisation'+str(runNumber)+'.root')
    trfCmplxFile = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-TrfLineCmplx-Params-v1.0.csv')
    print ("numSIMPATH, filename, rootinputfilename, rootfilename, trfCmplxFile \n", nuSIMPATH, "\n", filename, "\n", rootInputFilename, "\n", rootFilename,
         "\n", trfCmplxFile)
    outFilename = rootFilename
    logging.info("Parameters: %s,  \ntransfer line parameters: %s,  \noutput file: %s", filename,  trfCmplxFile, rootFilename)


# Get machine and run parameters
    RndmGen = Rndm.RandomGenerator(rootInputFilename,'histP','histXPS')
    nuStrt = nuPrdStrt.nuSTORMPrdStrght(filename)
    psLength = nuStrt.ProdStrghtLen()
    nuTrLnCmplx = nuTrf.nuSTORMTrfLineCmplx(trfCmplxFile)
    tlCmplxLength = nuTrLnCmplx.TrfLineCmplxLen()
    detectorPosZ = nuStrt.HallWallDist()
# set up the detector front face
    fluxPlane = plane.plane(psLength, detectorPosZ)
# set up the event history - instantiate
    eH = eventHistory.eventHistory()
    eH.outFile(outFilename)
# create the root structure to write to
    eH.rootStructure()
# initialise the python arrays to something sensible
    eH.makeHistory()

# event loop
for event in range(nEvents):
# set its values
    s = 0.0
    xl, xpl = RndmGen.getRandom2D()
    yl, ypl = RndmGen.getRandom2D()
    zl = -tlCmplxLength
    pPion = RndmGen.getRandom()
    pxl = pPion*xpl
    pyl = pPion*ypl
    pzl = np.sqrt(pPion*pPion - pxl**2 - pyl**2)
    t = nuTrLnCmplx.GenerateTime()*1E9
# transform to the global system
    xg, yg, zg, pxg, pyg, pzg = normInst.tltoGlbl(xl, yl, zl, pxl, pyl, pzl)
#  pion at target is a point source but with a momentum spread
    pionTarget = particle.particle(runNumber, event, s, xg, yg, zg, pxg, pyg, pzg, t, eventWeight, "pi+")
    eH.addParticle('target', pionTarget)

# generate a pion
    pi = piEvtInst.PionEventInstance(particleTar=pionTarget)
    tsc = pi.getTraceSpaceCoord()

# get the decay length - in the transfer line, in the production straight - or lost beyond the straight
    lifetime = pi.getLifetime()
    pathLength = lifetime*pPion*c/piMass

# decay in the transferline
    if ((tlFlag) and (pathLength < tlCmplxLength)):
      normInst.tlDecay()
    else:
# pion reaches end of transfer line just write out a new pion with altered s and z - all other pions must do this
# the magnets bend the beam so that the local co-ordinates are now the global ones
      se = tlCmplxLength
      ze = 0.0
# t = d/(beta*c)
      te = t + 1E9*se*math.sqrt(pPion**2 + piMass**2)/(c*pPion)
# x local is the same as before.
      pionPS = particle.particle(runNumber, event, se, xl, yl, ze, pxl, pyl, pzl, te, eventWeight, "pi+")
      eH.addParticle("productionStraight", pionPS)
# decay beyond the end of the production straight
    if ((lstFlag) and (pathLength > tlCmplxLength + psLength)):
        normInst.beyondPS()

# decay in the production straight
    if ((psFlag) and (pathLength >= tlCmplxLength) and (pathLength <= tlCmplxLength + psLength)):
        normInst.decayPiInPS()
# decay the muons
        if (muDcyFlag): normInst.decayMuons()


#  write to the root structure
    eH.fill()
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
    else:
        if (event%1000 ==0): print ("event number is ", event)


# Write to the root output file and close
eH.cd()
eH.write()
eH.outFileClose()
# Write out histograms
hm.histOutRoot("normalPlots.root")
##! Complete:
print()
print("========  Normalisation run : complete  ========")
