#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for muon class
=================================

  Assumes that nuSim code is in python path.

  Script tests the instantiation and get methods of the class

"""

import sys
import numpy as np
import math as mt
import muon as muon
import traceSpace
import MuonConst as MuonConst

##! Start:
classToTest="muon"
print("========  ", classToTest, ": tests start  ========")

##! Create instance, test built-in methods:
muonTest = 1
print()
print("muonTest:", muonTest, " Create muon and print quantities.")

mu = muon.muon(26, 134, 0.1, 0.15, -100.0, 0.0, 0.2, 0.22, 5.12, 0.00, 0.15)
print("    __str__:", mu)
print("    --repr__", repr(mu))
del mu

print()
muonTest = muonTest + 1
print("muonTest:", muonTest, " Create muon and check get methods.")
##! Create instance, test get methods
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

muCnst = MuonConst.MuonConst()
mass = muCnst.mass()/1000.

muErr = 0
mu = muon.muon(runNum, eventNum, s, x, y, z, px, py, pz, t, weight)

if mu.run() != runNum:
    muErr = muErr + 1
    print ("mu.run() is ", mu.run(), "and should be ", runNum)

if mu.event() != eventNum:
    muErr = muErr + 1
    print ("mu.event() is ", mu.event(), "and should be ", eventNum)

if mu.t() != t:
    muErr = muErr + 1
    print ("mu.t() is ", mu.t(), "and should be ", t)

if mu.weight() != weight:
    muErr = muErr + 1
    print ("mu.weight() is ", mu.weight(), "and should be ", weight)

if mu.mass() != mass:
    muErr = muErr + 1
    print ("mu.mass() is ", mu.mass(), "and should be ", mass)

if mu.traceSpace() != tSC:
    muErr = muErr + 1
    print ("mu.traceSpace() is ", mu.traceSpace(), "and should be ", tSC)

if mu.traceSpace().s() != s:
    muErr = muErr + 1
    print ("mu.traceSpace().s() is ", mu.traceSpace().s(), "and should be ", s)

if mu.traceSpace().x() != x:
    muErr = muErr + 1
    print ("mu.traceSpace().x() is ", mu.traceSpace().x(), "and should be ", x)

if mu.traceSpace().y() != y:
    muErr = muErr + 1
    print ("mu.traceSpace().y() is ", mu.traceSpace().y(), "and should be ", y)

if mu.traceSpace().z() != z:
    muErr = muErr + 1
    print ("mu.traceSpace().z() is ", mu.traceSpace().z(), "and should be ", z)

if mu.traceSpace().xp() != px/pz:
    muErr = muErr + 1
    print ("mu.traceSpace().xp() is ", mu.traceSpace().xp(), "and should be ", px/pz)

if mu.traceSpace().yp() != py/pz:
    muErr = muErr + 1
    print ("mu.traceSpace().yp() is ", mu.traceSpace().yp(), "and should be ", py/pz)

E = mt.sqrt(px*px + py*py + pz*pz + mass*mass)
pp = mu.p()
if pp[1][0] != px or pp[1][1] != py or pp[1][2] != pz or pp[0] != E:
    muErr = muErr + 1
    print ("mu.p() is ", mu.p(), "and should be ", fourVector)


##! Complete:
print()
print("========  muon:tests complete  ========")
print ("number of Errors is ", muErr)
sys.exit(0)