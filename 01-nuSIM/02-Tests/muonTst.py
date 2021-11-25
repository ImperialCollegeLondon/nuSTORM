#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for muon class
=================================

  Assumes that nuSim code is in python path.

  Script tests the instantiation and get methods of the class

"""

import sys, os
import numpy as np
import math as mt
import muon as muon
import traceSpace
import MuonConst as MuonConst

##! Start:
classToTest="muon"
print("========  ", classToTest, ": tests start  ========")
testFails = 0

##! Create instance, test built-in methods:
nuSIMPATH = os.getenv('nuSIMPATH')

##! Create instance, test built-in methods:
muonTest = 1
print()
print("muonTest:", muonTest, " Create muon and print quantities.")

mu = muon.muon(26, 134, 0.1, 0.15, -100.0, 0.0, 0.2, 0.22, 5.12, 0.00, 0.15)

restoreOut = sys.stdout
fileName = os.path.join(nuSIMPATH,'Scratch/muonTst.out')
refFile = os.path.join(nuSIMPATH,'02-Tests/referenceOutput/muonTst.ref')
sys.stdout = open(fileName,"w")
print("    __str__:", mu)
print("    --repr__", repr(mu))
sys.stdout.close()
sys.stdout = restoreOut
# compare the standard output to a reference file
a = os.popen('diff ' + fileName + ' ' + refFile)
output = a.read()
if (output == ""):
    pass
else:
    testFails = testFails + 1
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

mu = muon.muon(runNum, eventNum, s, x, y, z, px, py, pz, t, weight)

if mu.run() != runNum:
    testFails = testFails + 1
    print ("mu.run() is ", mu.run(), "and should be ", runNum)

if mu.event() != eventNum:
    testFails = testFails + 1
    print ("mu.event() is ", mu.event(), "and should be ", eventNum)

if mu.t() != t:
    testFails = testFails + 1
    print ("mu.t() is ", mu.t(), "and should be ", t)

if mu.weight() != weight:
    testFails = testFails + 1
    print ("mu.weight() is ", mu.weight(), "and should be ", weight)

if mu.mass() != mass:
    testFails = testFails + 1
    print ("mu.mass() is ", mu.mass(), "and should be ", mass)

if mu.traceSpace() != tSC:
    testFails = testFails + 1
    print ("mu.traceSpace() is ", mu.traceSpace(), "and should be ", tSC)

if mu.traceSpace().s() != s:
    testFails = testFails + 1
    print ("mu.traceSpace().s() is ", mu.traceSpace().s(), "and should be ", s)

if mu.traceSpace().x() != x:
    testFails = testFails + 1
    print ("mu.traceSpace().x() is ", mu.traceSpace().x(), "and should be ", x)

if mu.traceSpace().y() != y:
    testFails = testFails + 1
    print ("mu.traceSpace().y() is ", mu.traceSpace().y(), "and should be ", y)

if mu.traceSpace().z() != z:
    testFails = testFails + 1
    print ("mu.traceSpace().z() is ", mu.traceSpace().z(), "and should be ", z)

if mu.traceSpace().xp() != px/pz:
    testFails = testFails + 1
    print ("mu.traceSpace().xp() is ", mu.traceSpace().xp(), "and should be ", px/pz)

if mu.traceSpace().yp() != py/pz:
    testFails = testFails + 1
    print ("mu.traceSpace().yp() is ", mu.traceSpace().yp(), "and should be ", py/pz)

E = mt.sqrt(px*px + py*py + pz*pz + mass*mass)
pp = mu.p()
if pp[1][0] != px or pp[1][1] != py or pp[1][2] != pz or pp[0] != E:
    testFails = testFails + 1
    print ("mu.p() is ", mu.p(), "and should be ", fourVector)


##! Complete:
print()
print("========  muon:tests complete  ========")
print ("number of Errors is ", testFails)
if (testFails == 0):
    sys.exit(0)
else:
    sys.exit(1)