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
print("========  particle: tests start  ========")

##! Create instance, test built-in methods:
particleTest = 1
print()
print("particleTest:", particleTest, " Create particle and print quantities.")

p = particle.particle(0.1, 0.15, -100.0, 0.0, 0.2, 0.22, 5.12, 0.00, 0.16, 0.99)
print("    __str__:", p)
print("    --repr__", repr(p))
del p


print()
print("particleTest:", particleTest, " Create particle and check get methods.")
##! Create instance, test get methods
particleTest = particleTest + 1
print()
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
tSC = traceSpace.traceSpace(s, x, y, z, px/pz, py/pz)

parErr = 0
p = particle.particle(s, x, y, z, px, py, pz, t, weight, mass)

if p.t() != t:
    parErr = parErr + 1
    print ("p.t() is ", p.t(), "and should be ", t)

if p.weight() != weight:
    parErr = parErr + 1
    print ("p.weight() is ", p.weight(), "and should be ", weight)

if p.mass() != mass:
    parErr = parErr + 1
    print ("p.mass() is ", p.mass(), "and should be ", mass)

if p.traceSpace() != tSC:
    parErr = parErr + 1
    print ("p.traceSpace() is ", p.traceSpace(), "and should be ", tSC)

if p.traceSpace().s() != s:
    parErr = parErr + 1
    print ("p.traceSpace().s() is ", p.traceSpace().s(), "and should be ", s)

if p.traceSpace().x() != x:
    parErr = parErr + 1
    print ("p.traceSpace().x() is ", p.traceSpace().x(), "and should be ", x)

if p.traceSpace().y() != y:
    parErr = parErr + 1
    print ("p.traceSpace().y() is ", p.traceSpace().y(), "and should be ", y)

if p.traceSpace().z() != z:
    parErr = parErr + 1
    print ("p.traceSpace().z() is ", p.traceSpace().z(), "and should be ", z)

if p.traceSpace().xp() != px/pz:
    parErr = parErr + 1
    print ("p.traceSpace().xp() is ", p.traceSpace().xp(), "and should be ", px/pz)

if p.traceSpace().yp() != py/pz:
    parErr = parErr + 1
    print ("p.traceSpace().yp() is ", p.traceSpace().yp(), "and should be ", py/pz)

E = mt.sqrt(px*px + py*py + pz*pz + mass*mass)
pp = p.p()
if pp[1][0] != px or pp[1][1] != py or pp[1][2] != pz or pp[0] != E:
    parErr = parErr + 1
    print ("p.p() is ", p.p(), "and should be ", E, " ", px, "  ", py, " ", pz)


##! Complete:
print()
print("========  particle:tests complete  ========")
print ("number of Errors is ", parErr)
sys.exit(0)