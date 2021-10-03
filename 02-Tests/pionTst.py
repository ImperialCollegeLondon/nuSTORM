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
import pion as pion
import PionConst as PionConst
import traceSpace

##! Start:
print("========  pion: tests start  ========")

##! Create instance, test built-in methods:
pionTest = 1
print()
print("pionTest:", pionTest, " Create pion and print quantities.")

pi = pion.pion(26, 134, 0.1, 0.15, -100.0, 0.0, 0.2, 0.22, 5.12, 0.00, 0.16)
print("    __str__:", pi)
print("    --repr__", repr(pi))
del pi


print()
print("pionTest:", pionTest, " Create pion and check get methods.")
##! Create instance, test get methods
pionTest = pionTest + 1
print()
runNum = 137
eventNum = 695
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

piCnst = PionConst.PionConst()
mass = piCnst.mass()/1000.


piErr = 0
pi = pion.pion(runNum, eventNum, s, x, y, z, px, py, pz, t, weight)

if pi.run() != runNum:
    piErr = piErr + 1
    print ("pi.run() is ", pi.run(), "and should be ", runNum)

if pi.event() != eventNum:
    piErr = piErr + 1
    print ("pi.event() is ", pi.event(), "and should be ", eventNum)

if pi.t() != t:
    piErr = piErr + 1
    print ("pi.t() is ", pi.t(), "and should be ", t)

if pi.weight() != weight:
    piErr = piErr + 1
    print ("pi.weight() is ", pi.weight(), "and should be ", weight)

if pi.mass() != mass:
    piErr = piErr + 1
    print ("pi.mass() is ", pi.mass(), "and should be ", mass)

if pi.traceSpace() != tSC:
    piErr = piErr + 1
    print ("pi.traceSpace() is ", pi.traceSpace(), "and should be ", tSC)

if pi.traceSpace().s() != s:
    piErr = piErr + 1
    print ("pi.traceSpace().s() is ", pi.traceSpace().s(), "and should be ", s)

if pi.traceSpace().x() != x:
    piErr = piErr + 1
    print ("pi.traceSpace().x() is ", pi.traceSpace().x(), "and should be ", x)

if pi.traceSpace().y() != y:
    piErr = piErr + 1
    print ("pi.traceSpace().y() is ", pi.traceSpace().y(), "and should be ", y)

if pi.traceSpace().z() != z:
    piErr = piErr + 1
    print ("pi.traceSpace().z() is ", pi.traceSpace().z(), "and should be ", z)

if pi.traceSpace().xp() != px/pz:
    piErr = piErr + 1
    print ("pi.traceSpace().xp() is ", pi.traceSpace().xp(), "and should be ", px/pz)

if pi.traceSpace().yp() != py/pz:
    piErr = piErr + 1
    print ("pi.traceSpace().yp() is ", pi.traceSpace().yp(), "and should be ", py/pz)

E = mt.sqrt(px*px + py*py + pz*pz + mass*mass)
pp = pi.p()
if pp[1][0] != px or pp[1][1] != py or pp[1][2] != pz or pp[0] != E:
    piErr = piErr + 1
    print ("pi.p() is ", pi.p(), "and should be ", E, ", ", px, ", ", py, ", ", pz)


##! Complete:
print()
print("========  pion:tests complete  ========")
print ("number of Errors is ", piErr)
if (piErr == 0):
    sys.exit(0)
else:
    sys.exit(1)