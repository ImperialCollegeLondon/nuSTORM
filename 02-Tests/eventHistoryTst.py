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

##! addParticle location 'target' #############################################################################
descString = "Check addParticle and findParticle - with 'target' location"
descriptions.append(descString)

nTests = nTests + 1
print(testTitle, ": ",  descString)

obj = eventHistory.eventHistory()
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
testParticle = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, pdgCode)
obj.addParticle("target", testParticle)
par2 = obj.findParticle("target")
if (par2 == testParticle):
    pass
else:
    testFails = testFails + 1
    print(descriptions[nTests], " ..... failed")

del obj

##! addParticle location 'productionStraight' #############################################################################
descString = "Check addParticle and findParticle - with 'productionStraight' location"
descriptions.append(descString)

nTests = nTests + 1
print(testTitle, ": ",  descString)

obj = eventHistory.eventHistory()
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
testParticle = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, mass, pdgCode)
obj.addParticle("productionStraight", testParticle)
par2 = obj.findParticle("productionStraight")
if (par2 == testParticle):
    pass
else:
    testFails = testFails + 1
    print(descriptions[nTests], " ..... failed")

del obj


##! tests complete ########################################################################################

print()
print("========  ", testTitle, ": tests complete  ========")

print ("\nNumber of tests is ", nTests, " number of fails is ", testFails)

sys.exit(0)