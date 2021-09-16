#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for pion class
=================================

  Assumes that nuSim code is in python path.

  Script tests the instantiation and get methods of the class

"""

import sys
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

nTests = nTests + 1
print(testTitle, ": ",  descString)

p = particle.particle(42, 125, 0.1, 0.15, -100.0, 0.0, 0.2, 0.22, 5.12, 0.00, 0.16, 0.99, -14)
print("    __str__:", p)
print("    --repr__", repr(p))
del p


##! Create particle and check get methods #############################################################################
print()
descString = "Create particle check get methods"
descriptions.append(descString)

Tests = nTests + 1
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
mass = 0.99
pdgCode = 211
tSC = traceSpace.traceSpace(s, x, y, z, px/pz, py/pz)

p = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight, mass, pdgCode)

if p.run() != runNum:
    testFails = testFails + 1
    print ("p.run() is ", p.run(), "and should be ", runNum)

if p.event() != eventNum:
    testFails = testFails + 1
    print ("p.event() is ", p.event(), "and should be ", eventNum)

if p.t() != t:
    testFails = testFails + 1
    print ("p.t() is ", p.t(), "and should be ", t)

if p.weight() != weight:
    testFails = testFails + 1
    print ("p.weight() is ", p.weight(), "and should be ", weight)

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
    print ("p.traceSpace().xp() is ", p.traceSpace().xp(), "and should be ", px/pz)

if p.traceSpace().yp() != py/pz:
    testFails = testFails + 1
    print ("p.traceSpace().yp() is ", p.traceSpace().yp(), "and should be ", py/pz)

E = mt.sqrt(px*px + py*py + pz*pz + mass*mass)
pp = p.p()
if pp[1][0] != px or pp[1][1] != py or pp[1][2] != pz or pp[0] != E:
    testFails = testFails + 1
    print ("p.p() is ", p.p(), "and should be ", E, " ", px, "  ", py, " ", pz)


##! Test equality : ##################################################################
descString = "Checking equality for particles"
descriptions.append(descString)
print(testTitle, ": ",  descString)

p1 = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight, mass, pdgCode)

nTests=0
if (p == p1):
    pass
else:
    testFails = testFails + 1
    print(descriptions[nTests], " ..... failed")
    if p.run() != p1.run():
        print (".run() is ", p.run(), " and " , p1.run())
    if p.event() != p1.event():
        print ("p.event() is ", p.event(), " and " ,p1.event())
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
    if p.p() != p1.p():
        print ("p.p() is ", p.p(), " and " ,p1.p())

##! Test inequality: #################################################################
descString = "Checking inequality for particles"
descriptions.append(descString)
print(testTitle, ": ",  descString)

p2 = particle.particle(runNum, eventNum, s+1.0, x, y, z, px, py, pz, t, weight, mass, pdgCode)
nTests = nTests + 1
if (p != p2):
    pass
else:
    testFails = testFails + 1
    print(descriptions[nTests], " ..... failed")

##! Complete:
print()
print("========  particle:tests complete  ========")
print ("number of Errors is ", testFails)
sys.exit(0)