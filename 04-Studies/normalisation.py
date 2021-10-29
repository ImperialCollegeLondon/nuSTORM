#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model for calculating normalised numbers
========================================

    Assumes that nuSim code is in python path.

    @file normalisation.py

    @brief   calculates the normalisation of the data

    @author  Paul Kyberd

    @version     1.0
    @date        05 October 2021


"""

import os, sys
import numpy as np
import math as math
import PionEventInstance as piEvtInst
import NeutrinoEventInstance as nuEvtInst
import traceSpace as trSp
import particle as particle
import PionDecay as PionDecay
import PionConst as PC
import eventHistory as eventHistory

# calculate where in the ring, the decay occurred
def ring(pathLength):

  sectionLengths=[180.0, 150.0, 180.0, 150.0]
  arcRadius = 150.0/math.pi
#  print ("arcRadius is ", arcRadius)

  sectionPnt = 0
  while pathLength > 0.0:
    pathLength = pathLength - sectionLengths[sectionPnt]
    sectionPnt = (sectionPnt+1)%4
#    print ("pathlength is ", pathLength, "    sectionPnt is ", sectionPnt)

#  print ("pathlength is ", pathLength, "    sectionPnt is ", sectionPnt)
  sectionPnt = (sectionPnt - 1)%4
  decayPnt = pathLength + sectionLengths[sectionPnt]
#  print ("decay point is ", decayPnt, "    sectionPnt is ", sectionPnt)

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

printFlg = True
##! Start:
print("========  Narmalisation run: start  ========")
print()

#deltaX, zpos, direction = ring(179.9)
#print ("deltaX ", deltaX, "    zpos ", zpos, "   direction ", direction)
#deltaX, zpos, direction = ring(217.5)


# Get the nuSIM path name and use it to set names for the inputfile and the outputFile
nuSIMPATH = os.getenv('nuSIMPATH')
filename  = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
rootfilename = os.path.join(nuSIMPATH, 'Scratch/normalisation.root')

if (printFlg): print(".... file names set up")

#  Loop over events

piCnst  = PC.PionConst()
c = piCnst.SoL()
piMass = piCnst.mass()/1000.0
transferLine = 100.00                   # Transfer line length
runNumber = 50
nEvents = 6000
# set up the event history
eH = eventHistory.eventHistory()
eH.outFile("norm.root")
eH.rootStructure()
#  make sure the data structure is complete
testParticle = particle.particle(runNumber, -1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0,   "pi+")
eH.addParticle("target", testParticle)
eH.addParticle("productionStraight", testParticle)
eH.addParticle("pionDecay", testParticle)
eH.addParticle("muonProduction", testParticle)
eH.addParticle("piFlashNu", testParticle)
eH.addParticle("muonDecay", testParticle)
eH.addParticle("eProduction", testParticle)
eH.addParticle("numuProduction", testParticle)
eH.addParticle("nueProduction", testParticle)
eH.addParticle("numuDetector", testParticle)
eH.addParticle("nueDetector", testParticle)

crossSection = 50
eventWeight = crossSection/nEvents
#   Set the pion momentum
pionMom=6.0
nPass = 0
nDecay = 0


for event in range(nEvents):
# get a pion
    pi = piEvtInst.PionEventInstance(pionMom)
    if (event < 2): print ("pi is ", pi)
    tSC = pi.getTraceSpaceCoord()
    if (event < 2): print ("tSC is ", tSC)
    s = tSC[0]
    x = tSC[1]
    y = tSC[2]
    z = tSC[3]
    xp = tSC[4]
    yp = tSC[5]
    pPion = pi.getppiGen()
    px = pPion*xp
    py = pPion*yp
    pz = np.sqrt(pPion*pPion - px**2 - py**2)
    t = 0.0
#  pion at target has the same x,y, xp, yp as those at decay - but the s and z should be different 
    pionTarget = particle.particle(runNumber, event, 0.0, x, y, -100.0, px, py, pz, t, eventWeight, "pi+")
    eH.addParticle('target', pionTarget)
#

    if (event < 2): print ("pionProd is ", pionTarget)
# get check the decay length
    Dcy = PionDecay.PionDecay()
    lifetime = Dcy.getLifetime()
    pathLength = lifetime*pPion*c/piMass
    if (pathLength < transferLine):
# decays at present just set the productionStraight particle with a weight of zero - but updated s and z (small value of pz
# for tracespace calculation)
      noParticle = particle.particle(runNumber, event, transferLine, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, "pi+")
      eH.addParticle("productionStraight", noParticle)
# add a pion decay particle
      pionTLDecay = particle.particle(runNumber, event, s, x, y, z, px, py, pz, t, 0.0, "pi+")
      eH.addParticle("pionDecay", pionTLDecay)
      nDecay = nDecay + 1
    else:
# pion reaches end of transfer line just write out a new pion with altered s and z
      s = 100.0
      z = 0.0
      pionPS = particle.particle(runNumber, event, s, x, y, z, px, py, pz, t, eventWeight, "pi+")
      eH.addParticle("productionStraight", pionPS)
      nPass = nPass + 1

# now to find where in the ring the decay point occurred.
      deltaX, zpos, direction = ring(pathLength)
      s = 100.0+pathLength
      z = zpos
      x = x + deltaX
      pionDecay = particle.particle(runNumber, event, s, x, y, z, px, py, pz, t, eventWeight, "pi+")
      eH.addParticle("pionDecay", pionDecay)
# add the muon
      mu = pi.getmu4mmtm()
      if (event < 2): print("mu is ", mu)
      px = mu[1][0]
      py = mu[1][1]
      py = mu[1][2]
      muProd = particle.particle(runNumber, event, s, x, y, z, px, py, pz, t, eventWeight, "mu+")
      eH.addParticle("muonProduction", muProd)
# add the pion flash neutrino
      numu = pi.getnumu4mmtm()
      px = numu[1][0]
      py = numu[1][1]
      py = numu[1][2]
      if (event < 2): print("numu is ", numu)
      nuFlash = particle.particle(runNumber, event, s, x, y, z, px, py, pz, t, eventWeight, "numu")
      eH.addParticle("piFlashNu", nuFlash)
# decay the muon
      pbeam = 5.0
      nuEvt = nuEvtInst.NeutrinoEventInstance(pbeam)


# now find the decay point in the ring
    eH.fill()
# Make sure there is a structure - having zeroed the history after filling
    eH.addParticle("target", testParticle)
    eH.addParticle("productionStraight", testParticle)
    eH.addParticle("pionDecay", testParticle)
    eH.addParticle("muonProduction", testParticle)
    eH.addParticle("piFlashNu", testParticle)
    eH.addParticle("muonDecay", testParticle)
    eH.addParticle("eProduction", testParticle)
    eH.addParticle("numuProduction", testParticle)
    eH.addParticle("nueProduction", testParticle)
    eH.addParticle("numuDetector", testParticle)
    eH.addParticle("nueDetector", testParticle)


print ("pions which exit the transfer line ", nPass, "   those that decay ", nDecay, "  decay percent is ", nDecay/nEvents)
eH.write()
eH.outFileClose()

##! Complete:
print()
print("========  Normalisation run : complete  ========")
