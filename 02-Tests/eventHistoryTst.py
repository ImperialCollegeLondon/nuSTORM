#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for eventHistory class
=================================+

  Assumes that nuSim code is in python path.

    @file eventHistoryTst.py

    @brief   runBeamData programme
    Run and test the eventHistory class

    @author  Paul Kyberd

    @version     1.0
    @date        16 September 2021


"""

# Generic Python imports
import sys

from pathlib import Path            # checking for file existance

# nuStorm imports
import eventHistory as eventHistory
import particle as particle

##! Start:
nTests = 0
testFails = 0
descriptions=[]
locations=[]
testStatus=[]
testTitle = "eventHistory"

print("========  ", testTitle, ": tests start  ========")

##! Create instance, test built-in methods: #######################################################
descString = "Create instance and test built in methods"
descriptions.append(descString)

print(testTitle, ": ",  descString)

obj = eventHistory.eventHistory()
print("    __str__:", obj)
print("    --repr__", repr(obj))
del obj

##! outFilename check #############################################################################
descString = "Check output file is created"
descriptions.append(descString)

nTests = nTests + 1
print(testTitle, ": ",  descString)

obj = eventHistory.eventHistory()
obj.outFile("testFile.root")
obj.outFileClose()
file = Path("testFile.root")
if file.is_file():
    pass
else:
    testFails = testFails + 1
    print(descriptions[nTests], " ..... failed")
del obj

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

print (locations)
print (descriptions)

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
    testParticle = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, pdgCode)
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

##! display the event History #############################################################################
descString = "display the eventHistory"
descriptions.append(descString)

nTests = nTests + 1
print(testTitle, ": ",  descString)

obj = eventHistory.eventHistory()
testParticle = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, pdgCode)
testParticle1 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, -pdgCode)
testParticle2 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, pdgCode+10)
testParticle3 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, pdgCode+15)
testParticle4 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, 12)
testParticle5 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, 21)
testParticle6 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, -21)
testParticle7 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, -12)
testParticle8 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, 14)
testParticle9 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, 16)
testParticle10 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, -15)
obj.addParticle("target", testParticle)
obj.addParticle("productionStraight", testParticle1)
obj.addParticle("pionDecay", testParticle2)
obj.addParticle("muonProduction", testParticle3)
obj.addParticle("piFlashNu", testParticle4)
obj.addParticle("muonDecay", testParticle5)
obj.addParticle("eProduction", testParticle6)
obj.addParticle("numuProduction", testParticle7)
obj.addParticle("nueProduction", testParticle8)
obj.addParticle("numuDetector", testParticle9)
obj.addParticle("nueDetector", testParticle10)

obj.display()

del obj

##! tests complete ########################################################################################

print()
print("========  ", testTitle, ": tests complete  ========")

print ("\nNumber of tests is ", nTests, " number of fails is ", testFails)

sys.exit(0)