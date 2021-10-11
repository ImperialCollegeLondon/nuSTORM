#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for pion class
=================================

  Assumes that nuSim code is in python path.

  Script tests the instantiation and get methods of the class

Version history:
----------------------------------------------
 1.1: 21Sep21: Test for the new constructor
 1.0: 29Aug21: Test the class


"""

import sys
import os
import numpy as np
import math as mt
import particle as particle
import traceSpace

##! Start:

nTests = 0
testFails = 0
descriptions=[]
testStatus=[]
testTitle = "particle"

print("========  ", testTitle, ": tests start  ========")


##! Create instance, test built-in methods:
##! Create particle and print out #############################################################################
descString = "Create particle and print quantities"
descriptions.append(descString)

print(testTitle, ": ",  descString)

p = particle.particle(42, 125, 132.1, 0.1, 0.15, -100.0, 0.0, 0.2, 0.22, 5.12, 0.001, "pi+")
#    def __init__(self, runNum, eventNum, s, x, y, z, px, py, pz, t, eventWeight, particleType):
#  redirect standard out to a file
restoreOut = sys.stdout
sys.stdout = open("Scratch/particleTst.out","w")
print("    __str__:", p)
print("    --repr__", repr(p))
sys.stdout.close()
sys.stdout = restoreOut
# compare the standard output to a reference file
a = os.popen('diff Scratch/particleTst.out 02-Tests/particleTst.ref')
output = a.read()
if (output == ""):
    pass
else:
    testFails = testFails + 1
del p

##! Create particle and check get methods #############################################################################
nTests = nTests + 1
descString = "Create particle check get methods"
descriptions.append(descString)

print(testTitle, ": ",  descString)

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
weight = 0.15
tSC = traceSpace.traceSpace(s, x, y, z, px/pz, py/pz)

p = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight, "pi+")

if p.run() != runNum:
    testFails = testFails + 1
    print ("p.run() is ", p.run(), "and should be ", runNum)

if p.event() != eventNum:
    testFails = testFails + 1
    print ("p.event() is ", p.event(), "and should be ", eventNum)

pdgCode = 211
if p.pdgCode() != pdgCode:
    testFails = testFails + 1
    print ("p.pdgCode() is ", p.pdgCode(), "and should be ", pdgCode)

if p.t() != t:
    testFails = testFails + 1
    print ("p.t() is ", p.t(), "and should be ", t)

if p.weight() != weight:
    testFails = testFails + 1
    print ("p.weight() is ", p.weight(), "and should be ", weight)

mass = 139.57061/1000.0
if p.mass() != mass:
    testFails = testFails + 1
    print ("p.mass() is ", p.mass(), "and should be ", mass)

if p.traceSpace() != tSC:
    testFails = testFails + 1
    print ("p.traceSpace() is ", p.traceSpace(), "and should be ", tSC)

if p.traceSpace().s() != s:
    testFails = testFails + 1
    print ("p.traceSpace().s() is ", p.traceSpace().s(), "and should be ", s)

if p.traceSpace().x() != x:
    testFails = testFails + 1
    print ("p.traceSpace().x() is ", p.traceSpace().x(), "and should be ", x)

if p.traceSpace().y() != y:
    testFails = testFails + 1
    print ("p.traceSpace().y() is ", p.traceSpace().y(), "and should be ", y)

if p.traceSpace().z() != z:
    testFails = testFails + 1
    print ("p.traceSpace().z() is ", p.traceSpace().z(), "and should be ", z)

if p.traceSpace().xp() != px/pz:
    testFails = testFails + 1
    print ("p.traceSpace().xp() is ", p.traceSpace().yp(), "and should be ", px/pz)

if p.traceSpace().yp() != py/pz:
    testFails = testFails + 1
    print ("p.traceSpace().yp() is ", p.traceSpace().yp(), "and should be ", py/pz)

E = mt.sqrt(px*px + py*py + pz*pz + mass*mass)
pp = p.p()
if pp[1][0] != px or pp[1][1] != py or pp[1][2] != pz or pp[0] != E:
    testFails = testFails + 1
    print ("p.p() is ", p.p(), "and should be ", E, " ", px, "  ", py, " ", pz)


##! Test equality : ##################################################################
nTests = nTests + 1
descString = "Checking equality for particles"
descriptions.append(descString)
print(testTitle, ": ",  descString)

p1 = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight, "pi+")

if (p == p1):
    pass
else:
    testFails = testFails + 1
    print(descriptions[nTests], " ..... failed")
    if p.run() != p1.run():
        print (".run() is ", p.run(), " and " , p1.run())
    if p.event() != p1.event():
        print ("p.event() is ", p.event(), " and " ,p1.event())
    if p.pdgCode() != p1.pdgCode():
         print ("p.pdgCode() is ", p.pdgCode(), " and " ,p1.pdgCode())       
    if p.s() != p1.s():
        print ("p.s() is ", p.s(), " and " ,p1.s())
    if p.x() != p1.x():
        print ("p.x() is ", p.x(), " and " ,p1.x())
    if p.y() != p1.y():
        print ("p.y() is ", p.y(), " and " ,p1.y())
    if p.z() != p1.z():
        print ("p.z() is ", p.z(), " and " ,p1.z())
    if p.xp() != p1.xp():
        print ("p.xp() is ", p.xp(), " and " ,p1.xp())
    if p.yp() != p1.yp():
        print ("p.yp() is ", p.yp(), " and " ,p1.yp())
    if p.t() != p1.t():
        print ("p.t() is ", p.t(), " and " ,p1.t())
    if p.weight() != p1.weight():
        print ("p.weight() is ", p.weight(), " and " ,p1.weight())
    if p.mass() != p1.mass():
        print ("p.mass() is ", p.mass(), " and " ,p1.mass())
    if p.p()[0] != p1.p()[0]:
        print ("p.p()[0] is ", p.p()[0], " and " ,p1.p()[0])
    if p.p()[1][0] != p1.p()[1][0]:
        print ("p.p()[1][0] is ", p.p()[1][0], " and " ,p1.p()[1][0])
    if p.p()[1][1] != p1.p()[1][1]:
        print ("p.p()[1][1] is ", p.p()[1][1], " and " ,p1.p()[1][1])
    if p.p()[1][2] != p1.p()[1][2]:
        print ("p.p()[1][2] is ", p.p()[1][2], " and " ,p1.p()[1][2])

##! Test inequality: #################################################################
nTests = nTests + 1
descString = "Checking inequality for particles"
descriptions.append(descString)
print(testTitle, ": ",  descString)

tstFlag = True

pa = particle.particle(runNum+1, eventNum, s, x, y, z, px, py, pz, t, weight, "mu-")
if (p != pa):
    pass
else:
    tstFlag = False
pb = particle.particle(runNum, eventNum+1, s, x, y, z, px, py, pz, t, weight, "mu-")
if (p != pb):
    pass
else:
    tstFlag = False
pc = particle.particle(runNum, eventNum, s+1.0, x, y, z, px, py, pz, t, weight, "mu-")
if (p != pc):
    pass
else:
    tstFlag = False
pd = particle.particle(runNum, eventNum, s, x+1.0, y, z, px, py, pz, t, weight, "mu-")
if (p != pd):
    pass
else:
    tstFlag = False
pe = particle.particle(runNum, eventNum, s, x, y+1.0, z, px, py, pz, t, weight, "mu-")
if (p != pe):
    pass
else:
    tstFlag = False
pf = particle.particle(runNum, eventNum, s, x, y, z+1.0, px, py, pz, t, weight, "mu-")
if (p != pf):
    pass
else:
    tstFlag = False
pg = particle.particle(runNum, eventNum, s, x, y, z, px+1.0, py, pz, t, weight, "mu-")
if (p != pg):
    pass
else:
    tstFlag = False
ph = particle.particle(runNum, eventNum, s, x, y, z, px, py+1.0, pz, t, weight, "mu-")
if (p != ph):
    pass
else:
    tstFlag = False
pi = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz+1.0, t, weight, "mu-")
if (p != pi):
    pass
else:
    tstFlag = False
pj = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t+1.0, weight, "mu-")
if (p != pj):
    pass
else:
    tstFlag = False
pk = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight+1.0, "mu-")
if (p != pk):
    pass
else:
    tstFlag = False
pl = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight, "e-")
if (p != pl):
    pass
else:
    tstFlag = False


if (tstFlag == False):
    testFails = testFails + 1
    print(descriptions[nTests], " ..... failed")

##! Test new constructor: #################################################################
nTests = nTests + 1
descString = "Checking new constructor"
descriptions.append(descString)
print(testTitle, ": ",  descString)

px = 1.02
py = 2.04
pz = 9.18
part = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight, "e+")
mass = part.mass()
E = mt.sqrt(mass*mass + px*px + py*py + pz*pz)
p3 = np.array([px, py, pz])
p4 = np.array([E, p3],dtype=object)
altPart = particle.particle(runNum, eventNum, s, x, y, z, p4, t, weight, "e+")

if (part == altPart):
    pass
else:
    testFails = testFails + 1
    print(descriptions[nTests], " ..... failed")

##! Complete:
print()
print("========  particle:tests complete  ========")
print ("\nNumber of tests is ", nTests, " number of fails is ", testFails)
if testFails == 0:
    sys.exit(0)
else:
    sys.exit(1)