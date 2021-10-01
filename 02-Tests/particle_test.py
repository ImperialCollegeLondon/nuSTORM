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


def test_print():
    nuSIMPATH = os.getenv('nuSIMPATH')
##! Create instance, test built-in methods:
##! Create particle and print out #############################################################################
    descString = "Create particle and print quantities"
    descriptions.append(descString)

    print(testTitle, ": ",  descString)

    p = particle.particle(42, 125, 132.1, 0.1, 0.15, -100.0, 0.0, 0.2, 0.22, 5.12, 0.001, "pi+")
#    def __init__(self, runNum, eventNum, s, x, y, z, px, py, pz, t, eventWeight, particleType):
#  redirect standard out to a file
    restoreOut = sys.stdout
    fileName = os.path.join(nuSIMPATH,'Scratch/particleTst.out')
    refFile = os.path.join(nuSIMPATH,'02-Tests/referenceOutput/particleTst.ref')
    print ('filename is ', fileName)
    print ('refFile is ', refFile)
    sys.stdout = open(fileName,"w")
    print("    __str__:", p)
    print("    --repr__", repr(p))
    sys.stdout.close()
    sys.stdout = restoreOut
# compare the standard output to a reference file
    a = os.popen('diff ' + fileName + ' ' + refFile)
    output = a.read()
    assert output == ""
    if (output == ""):
        pass
    else:
        testFails = testFails + 1
    del p

##! Create particle and check get methods #############################################################################
def test_get():
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

    assert p.run() == runNum
    assert p.event() == eventNum
    pdgCode = 211
    assert p.pdgCode() == pdgCode
    assert p.t() == t
    assert p.weight() == weight
    mass = 139.57061/1000.0
    assert p.mass() == mass
    assert p.traceSpace() == tSC
    assert p.traceSpace().s() == s
    assert p.traceSpace().x() == x
    assert p.traceSpace().y() == y
    assert p.traceSpace().z() == z
    assert p.traceSpace().xp() == px/pz
    assert p.traceSpace().yp() == py/pz

    E = mt.sqrt(px*px + py*py + pz*pz + mass*mass)
    assert p.p()[1][0] == px
    assert p.p()[1][1] == py
    assert p.p()[1][2] == pz
    assert p.p()[0] == E

##! Test equality : ##################################################################
def test_equal():
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

    descString = "Checking equality for particles"
    descriptions.append(descString)
    print(testTitle, ": ",  descString)

    p = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight, "pi+")
    p1 = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight, "pi+")

    if (p == p1):
        pass
    else:
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
        assert p == p1

##! Test inequality : ##################################################################
def test_inequal():
    descString = "Checking inequality for particles"
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

    pa = particle.particle(runNum+1, eventNum, s, x, y, z, px, py, pz, t, weight, "pi+")
    assert p != pa
    pb = particle.particle(runNum, eventNum+1, s, x, y, z, px, py, pz, t, weight, "pi+")
    assert p != pb
    pc = particle.particle(runNum, eventNum, s+1.0, x, y, z, px, py, pz, t, weight, "pi+")
    assert p != pc
    pd = particle.particle(runNum, eventNum, s, x+1.0, y, z, px, py, pz, t, weight, "pi+")
    assert p != pd
    pe = particle.particle(runNum, eventNum, s, x, y+1.0, z, px, py, pz, t, weight, "pi+")
    assert p != pe
    pf = particle.particle(runNum, eventNum, s, x, y, z+1.0, px, py, pz, t, weight, "pi+")
    assert p != pf
    pg = particle.particle(runNum, eventNum, s, x, y, z, px+1.0, py, pz, t, weight, "pi+")
    assert p != pg
    ph = particle.particle(runNum, eventNum, s, x, y, z, px, py+1.0, pz, t, weight, "pi+")
    assert p != ph
    pi = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz+1.0, t, weight, "pi+")
    assert p != pi
    pj = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t+1.0, weight, "pi+")
    assert p != pj
    pk = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight+1.0, "pi+")
    assert p != pk
    pl = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight, "e-")
    assert p != pl

##! Test new constructor: #################################################################
def test_construct():

    descString = "Checking new constructor"
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
    px = 1.02
    py = 2.04
    pz = 9.18
    part = particle.particle(runNum, eventNum, s, x, y, z, px, py, pz, t, weight, "e+")
    mass = part.mass()
    E = mt.sqrt(mass*mass + px*px + py*py + pz*pz)
    p3 = np.array([px, py, pz])
    p4 = np.array([E, p3],dtype=object)
    altPart = particle.particle(runNum, eventNum, s, x, y, z, p4, t, weight, "e+")
    assert part == altPart

##! Complete:
print()
print("========  particle:tests complete  ========")
print ("\nNumber of tests is ", nTests, " number of fails is ", testFails)
