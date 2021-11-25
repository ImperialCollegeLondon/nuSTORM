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
import csv                          # so I can read a synthetic data file

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

##! rRead and write the event history tree: #######################################################
descString = "Read and Write eventHistory tree"
descriptions.append(descString)

print(testTitle, ": ",  descString)

# open the root output file
obj = eventHistory.eventHistory()
obj.outFile("testFile.root")
obj.rootStructure()

# open the file of dummy data
s = 10.8
px = 0.82
py = 0.45
pz = 4.67
t = 24.5
eventWeight = 0.15
mass = 0.99
pdgCode = 211

dummyData = open('02-Tests/dummyData.csv')
readDummy = csv.reader(dummyData)
header=[]
header=next(readDummy)

for instance in readDummy:
    loc = instance[0]
    if loc == "endEvent":
        obj.fill()
    else:
        runNum = int(instance[1])
        eventNum = int(instance[2])
        x = float(instance[3])
        y = float(instance[4])
        z = float(instance[5])
        s = float(instance[6])
        px = float(instance[7])
        py = float(instance[8])
        pz = float(instance[9])
        t = float(instance[10])
        eventWeight = float(instance[11])
        pType = instance[12]

        testParticle = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, eventWeight, pType)
        obj.addParticle(loc, testParticle)
#        obj.display()

# write out data and lose the outfile
obj.write()
obj.outFileClose()

print("finished writing")
objRd = eventHistory.eventHistory()
objRd.inFile("testFile.root")

nEvent = objRd.getEntries()

#   rewind the input csv file and read it again
dummyData.seek(0)
#   first the header
header=next(readDummy)

testFlag = True

print("start reading")
# Now read the input from the csv file, create the eventHistory and com[are it to the one recovered from the root file
#  get the first event fromnext event from root
locRecord=[]
locStatus=[]

objRd.readNext()
for instance in readDummy:
    loc = instance[0]
    if loc == "endEvent":
        objRd.readNext()
    else:
        runNum = int(instance[1])
        eventNum = int(instance[2])
        x = float(instance[3])
        y = float(instance[4])
        z = float(instance[5])
        s = float(instance[6])
        px = float(instance[7])
        py = float(instance[8])
        pz = float(instance[9])
        t = float(instance[10])
        eventWeight = float(instance[11])
        pType = instance[12]

#   we have the information for a particle at a location from the csv file use it to create a particle
        csvPart = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, eventWeight, pType)
#   compare that particle with the particle from the corresponding location in the root event
        rootPart = objRd.findParticle(loc)
        locRecord.append(loc)
        if (csvPart == rootPart):
            locStatus.append(True)
        else:
            testFlag = False
            locStatus.append(False)

if testFlag == False:
    print("testflag is ", testFlag)
    testFails = testFails + 1

for pnt in range(11):
    print("status is ", locStatus[pnt], "     position is ", locRecord[pnt])


del obj
##! Create instance, test built-in methods: #######################################################
descString = "Create instance and test built in methods"
descriptions.append(descString)

nTests = nTests + 1
print(testTitle, ": ",  descString)

obj = eventHistory.eventHistory()
print("    __str__:", obj)
print("    --repr__", repr(obj))
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

##! display the event History #############################################################################
descString = "display the eventHistory"
descriptions.append(descString)

nTests = nTests + 1
print(testTitle, ": ",  descString)

obj = eventHistory.eventHistory()
testParticle = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight,   "pi+")
testParticle1 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight,  "pi-")
testParticle2 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight,  "mu+")
testParticle3 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight,  "mu-")
testParticle4 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight,  "e+")
testParticle5 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight,  "e-")
testParticle6 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight,  "nue")
testParticle7 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight,  "nueBar")
testParticle8 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight,  "numu")
testParticle9 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight,  "numuBar")
testParticle10 = particle.particle(runNum, eventNum, x, y, z, s, px, py, pz, t, eventWeight, "pi+")
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