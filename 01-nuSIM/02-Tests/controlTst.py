#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for control class
=============================

  Assumes that nuSim code is in python path.


Version history:
----------------------------------------------
 1.0: 07Jan22: Test the class

 1.1: 02Dec:    Add test for flukaFile directory


"""

import sys
import os
import numpy as np
import control

##! Start:

nTests = 0
testFails = 0
descriptions=[]
testStatus=[]
testTitle = "control"

print(f"========  {testTitle} : tests start  ========")


##! Create instance, test built-in methods:
##! Create particle and print out #############################################################################
nTests = nTests + 1
descString = "Create control - and dump contents of controlFile"
descriptions.append(descString)
print(f"{testTitle} :  {descString}")

con = control.control("02-Tests/referenceOutput/PiFlash40-304.dict")

print(f"con is {con} ... ")
try:
    print(f"flukaFileDir is {con.flukaFileDir()} ... ")
except KeyError:
    print("key error thrown")

##! Create particle and check get methods #############################################################################
nTests = nTests + 1
descString = "Check normalisation process flags"
descriptions.append(descString)
if (con.tlFlag()):
    testStatus.append("True")
else:
    testFails = testFails + 1
if (con.psFlag()):
    print (f"psFlag is true")
else:
    print (f"psFlag is false")
if (con.lstFlag()):
    print (f"lstFlag is true")
else:
    print (f"lstFlag is false")
if (con.muDcyFlag()):
    print (f"muDcyFlag is true")
else:
    print (f"muDcyFlag is false")


##! Test equality : ##################################################################
nTests = nTests + 1
descString = "Checking run and event"
descriptions.append(descString)
print(f"{testTitle} :   {descString}")

print(f"runNumber is {con.runNumber()}")

print(f"runNumber is {con.runNumber(True)}")


print(f"number of Events to generate is {con.nEvents()}")


##! Test print flags : ###############################################################

print(f"\n\n========= Print flags check")

print(f"main routine print flag {con.mainPrnt()}")
print(f"neutrinoEventInstance routine print flag {con.nuEvtInstPrnt()}")


##! Complete:
print()
print(f"========  {testTitle}:tests complete  ========")
print (f"\nNumber of tests is {nTests} number of fails is {testFails}")
if testFails == 0:
    sys.exit(0)
else:
    sys.exit(1)