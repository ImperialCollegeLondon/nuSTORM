#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to analyset he eventHistory for a run
============================================

  Assumes that nuSim code is in python path.

    @file analyseEH.py

    @brief   runBeamData programme
    Run and test the eventHistory class

    @author  Paul Kyberd

    @version     1.0
    @date        29 October 2021


"""

# Generic Python imports
import sys

from pathlib import Path            # checking for file existance
import csv                          # so I can read a synthetic data file

# nuStorm imports
import eventHistory as eventHistory
import particle as particle
import histoManager
import histsCreate

##! Start:
nTests = 0
testFails = 0
descriptions=[]
locations=[]
testStatus=[]
testTitle = "analyse Event History"

print("========  ", testTitle, ": start  ========")

##! rRead and write the event history tree: #######################################################
descString = "Reading an eventHistory tree"
descriptions.append(descString)

print(testTitle, ": ",  descString)

# Create histogram manager
hm = histoManager.histoManager()
hC=histsCreate.histsCreate(hm)
# book histograms
hTitle = "x position at Detector" 
hBins  = 100
hLower = -10.0
hUpper = 10.0
h1 = hm.book(hTitle, hBins, hLower, hUpper)
#
hTitle = "s Prod Strght: no decay" 
hBins  = 100
hLower = -10.0
hUpper = 100.0
h2 = hm.book(hTitle, hBins, hLower, hUpper)
#
hTitle = "z Prod Strght: no decay" 
hBins  = 100
hLower = -10.0
hUpper = 100.0
h3 = hm.book(hTitle, hBins, hLower, hUpper)
#
hTitle = "t Prod Strght: no decay" 
hBins  = 100
hLower = 20.0
hUpper = 60.0
h4 = hm.book(hTitle, hBins, hLower, hUpper)
#
hTitle = "s Prod Strght: decay" 
hBins  = 100
hLower = -10.0
hUpper = 100.0
h5 = hm.book(hTitle, hBins, hLower, hUpper)
#
hTitle = "z Prod Strght: decay" 
hBins  = 100
hLower = -10.0
hUpper = 100.0
h6 = hm.book(hTitle, hBins, hLower, hUpper)
#
hTitle = "t Prod Strght: decay" 
hBins  = 100
hLower = 20.0
hUpper = 60.0
h7 = hm.book(hTitle, hBins, hLower, hUpper)
#
hTitle = "piPS weight" 
hBins  = 100
hLower = -1.0
hUpper = 99.0
h8 = hm.book(hTitle, hBins, hLower, hUpper)



#  create an eventHistory object
objRd = eventHistory.eventHistory()
#  set the file to read frp,
objRd.inFile("norm.root")

nEvent = objRd.getEntries()
print ("number of entries is ", nEvent)


hC.histAdd("target")
hC.histAdd("productionStraight")
hC.histAdd("pionDecay")
hC.histAdd("piFlashNu")
hC.histAdd("muonDecay")
hC.histAdd("eProduction")
hC.histAdd("numuProduction")
hC.histAdd("nueProduction")
hm.setParams("m")
testFlag = True

nTarget = 0
nProductionStraight = 0
nPionDecay = 0
for pnt in range(nEvent):
    objRd.readNext()
#  unpack required particles
    loc = "target"
    piTarget = objRd.findParticle(loc)
    hC.histsFill(loc, piTarget)
    if piTarget.event() > 0:
        nTarget = nTarget + 1
    loc = "productionStraight"
    piPS = objRd.findParticle(loc)
    hC.histsFill(loc, piPS)
    h8.Fill(piPS.weight())
    if piPS.weight() > 0:
        nProductionStraight = nProductionStraight + 1
        h2.Fill(piPS.s())
        h3.Fill(piPS.z())
        h4.Fill(piPS.t())
    else:
        h5.Fill(piPS.s())
        h6.Fill(piPS.z())
        h7.Fill(piPS.t())
    loc = "pionDecay"
    piDecay = objRd.findParticle(loc)
    hC.histsFill(loc, piDecay)
    if piPS.weight() > 0:
        nPionDecay = nPionDecay + 1
    loc = "piFlashNu"
    nuMuFlash = objRd.findParticle(loc)
    hC.histsFill(loc, nuMuFlash)




    loc = "numuDetector"
    numu = objRd.findParticle(loc)
#  unpack variables - only ploy if event number is non-negative
    if numu.event() > 0:
        x = numu.x()
        h1.Fill(x)

print ("number of pions at target ", nTarget, ":   number at production straight ", nProductionStraight, ":   number of pion decays ", nPionDecay)

# print out results
hm.histdo()


#eH.addParticle("target", testParticle)
#eH.addParticle("productionStraight", testParticle)
#eH.addParticle("pionDecay", testParticle)
#eH.addParticle("muonProduction", testParticle)
#eH.addParticle("piFlashNu", testParticle)
#eH.addParticle("muonDecay", testParticle)
#eH.addParticle("eProduction", testParticle)
#eH.addParticle("numuProduction", testParticle)
#eH.addParticle("nueProduction", testParticle)
##eH.addParticle("numuDetector", testParticle)
#eH.addParticle("nueDetector", testParticle)
#objRd.display();
#loc ="target"
#rootPart = objRd.findParticle(loc)
#print ("root part is ", rootPart)

'''
##! check addParticle with the different locations ########################################################################
locTest = 0
locations.append('target')
descString = "Check addParticle and findParticle - with '"+locations[-1]+" 'location"
descriptions.append(descString)
#
locations.append('productionStraight')
descString = "Check addParticle and findParticle - with '"+locations[-1]+"' location"
descriptions.append(descString)
#
locations.append('pionDecay')
descString = "Check addParticle and findParticle - with '"+locations[-1]+"' location"
descriptions.append(descString)
#
locations.append('muonProduction')
descString = "Check addParticle and findParticle - with '"+locations[-1]+"' location"
descriptions.append(descString)
#
locations.append('piFlashNu')
descString = "Check addParticle and findParticle - with '"+locations[-1]+"' location"
descriptions.append(descString)
#
locations.append('muonDecay')
descString = "Check addParticle and findParticle - with '"+locations[-1]+"' location"
descriptions.append(descString)
#
locations.append('eProduction')
descString = "Check addParticle and findParticle - with '"+locations[-1]+"' location"
descriptions.append(descString)
#
locations.append('numuProduction')
descString = "Check addParticle and findParticle - with '"+locations[-1]+"' location"
descriptions.append(descString)
#
locations.append('nueProduction')
descString = "Check addParticle and findParticle - with '"+locations[-1]+"' location"
descriptions.append(descString)
#
locations.append('numuDetector')
descString = "Check addParticle and findParticle - with '"+locations[-1]+"' location"
descriptions.append(descString)
#
locations.append('nueDetector')
descString = "Check addParticle and findParticle - with '"+locations[-1]+"' location"
descriptions.append(descString)

runNum = 43
eventNum = 67
x = 1.4
y = 2.3
z = 10.3
s = 10.8
px = 0.82
py = 0.45
pz = 4.67
t = 24.5
eventWeight = 0.15
mass = 0.99
pdgCode = 211

offset = 2
for pnt in range(len(locations)):
    nTests = nTests + 1
    print(testTitle, ": ",  descriptions[pnt+offset])
    obj = eventHistory.eventHistory()
    testParticle = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, pdgCode)
    obj.addParticle(locations[pnt], testParticle)
    retPar = obj.findParticle(locations[pnt])
    if (retPar == testParticle):
        pass
    else:
        testFails = testFails + 1
        print(descriptions[nTests+offset], " ..... failed")

del obj
del testParticle
del retPar
'''


print()
print("========  ", testTitle, ": complete  ========")

sys.exit(0)