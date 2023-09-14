#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for analyse class
=============================

  Assumes that nuSim code is in python path.


Version history:
----------------------------------------------
 1.0: 11Sep23: Test the class built in methods


"""

import sys
import os
import numpy as np
import analyse

##! Start:

nTests = 0
testFails = 0
descriptions=[]
testStatus=[]
testTitle = "analyse"

print(f"========  {testTitle} : tests start  ========")


##! Create instance, test built-in methods:
##! Create particle and print out #############################################################################
nTests = nTests + 1
descString = "Create analyse - and run built in methods - check by hand"
descriptions.append(descString)
print(f"{testTitle} :  {descString}")

an = analyse.analyse()

print(f"an is {an} ... ")

print(f"__repr__ is {repr(an)} ")

for pnt in range(100):
    for pnt1 in range(pnt):
        an.process(pnt)

an.conclude("analyseTest.root")

##! Complete:
print()
print(f"========  {testTitle}:tests complete  ========")
print (f"\nNumber of tests is {nTests} number of fails is {testFails}")
if testFails == 0:
    sys.exit(0)
else:
    sys.exit(1)