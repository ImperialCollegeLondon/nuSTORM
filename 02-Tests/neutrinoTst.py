#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for neutrino class
==============================

  Assumes that nuSim code is in python path.

  Script tests the instantiation and get methods of the class

"""

import sys
import numpy as np
import math as mt
import neutrino as neutrino
import traceSpace

##! Start:
print("========  neutrino: tests start  ========")

##! Create instance, test built-in methods:
neuTest = 1
print()
print("neutrino Test:", neuTest, " Create neutrino and print quantities.")

nu = neutrino.neutrino(26, 134, 0.1, 0.15, -100.0, 0.0, 0.2, 0.22, 5.12, 0.00, 0.16)
print("    __str__:", nu)
print("    --repr__", repr(nu))
del nu


print()
print("neutrino Test:", neuTest, " Create neutrino and check get methods.")
##! Create instance, test get methods
neuTest = neuTest + 1
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

mass = 0.0

neuErr = 0
neu = neutrino.neutrino(runNum, eventNum, s, x, y, z, px, py, pz, t, weight)

if neu.run() != runNum:
    neuErr = neuErr + 1
    print ("neu.run() is ", neu.run(), "and should be ", runNum)

if neu.event() != eventNum:
    neuErr = neuErr + 1
    print ("neu.event() is ", neu.event(), "and should be ", eventNum)

if neu.t() != t:
    neuErr = neuErr + 1
    print ("neu.t() is ", neu.t(), "and should be ", t)

if neu.weight() != weight:
    neuErr = neuErr + 1
    print ("neu.weight() is ", neu.weight(), "and should be ", weight)

if neu.mass() != mass:
    neuErr = neuErr + 1
    print ("neu.mass() is ", neu.mass(), "and should be ", mass)

if neu.traceSpace() != tSC:
    neuErr = neuErr + 1
    print ("neu.traceSpace() is ", neu.traceSpace(), "and should be ", tSC)

if neu.traceSpace().s() != s:
    neuErr = neuErr + 1
    print ("neu.traceSpace().s() is ", neu.traceSpace().s(), "and should be ", s)

if neu.traceSpace().x() != x:
    neuErr = neuErr + 1
    print ("neu.traceSpace().x() is ", neu.traceSpace().x(), "and should be ", x)

if neu.traceSpace().y() != y:
    neuErr = neuErr + 1
    print ("neu.traceSpace().y() is ", neu.traceSpace().y(), "and should be ", y)

if neu.traceSpace().z() != z:
    neuErr = neuErr + 1
    print ("neu.traceSpace().z() is ", neu.traceSpace().z(), "and should be ", z)

if neu.traceSpace().xp() != px/pz:
    neuErr = neuErr + 1
    print ("neu.traceSpace().xp() is ", neu.traceSpace().xp(), "and should be ", px/pz)

if neu.traceSpace().yp() != py/pz:
    neuErr = neuErr + 1
    print ("neu.traceSpace().yp() is ", neu.traceSpace().yp(), "and should be ", py/pz)

E = mt.sqrt(px*px + py*py + pz*pz + mass*mass)
pp = neu.p()
if pp[1][0] != px or pp[1][1] != py or pp[1][2] != pz or pp[0] != E:
    neuErr = neuErr + 1
    print ("neu.p() is ", neu.p(), "and should be ", fourVector)


##! Complete:
print()
print("========  neutrino :tests complete  ========")
print ("number of Errors is ", neuErr)
sys.exit(0)