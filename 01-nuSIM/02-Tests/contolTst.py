#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for control class
=============================

  Assumes that nuSim code is in python path.


Version history:
----------------------------------------------
 1.0: 07Jan22: Test the class


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

con = control.control()

print(f"con is {con} ... ")


##! Create particle and check get methods #############################################################################
nTests = nTests + 1
descString = "Check normalisation process flags"
descriptions.append(descString)
if (con.tlFlag()):
    print (f"tlFlag is true")
else:
    print (f"tlFlag is false")
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
print(f"number of Events to generate is {con.nEvents()}")


##! Complete:
print()
print(f"========  {testTitle}:tests complete  ========")
print (f"\nNumber of tests is {nTests} number of fails is {testFails}")
if testFails == 0:
    sys.exit(0)
else:
    sys.exit(1)