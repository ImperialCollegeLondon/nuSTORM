#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "PionConst" class ... initialisation and get methods
=================================

  PionConst.py -- set "relative" path to code

"""

import PionConst as pc
import sys, os

testFails = 0
##! Start:
print("========  PionConst: tests start  ========")


##! Test singleton class feature:
PionConstTest = 1
print()
print("PionConstTest:", PionConstTest, " check if class is a singleton.")
piCnst  = pc.PionConst()
piCnst1 = pc.PionConst()
print("    piCnst singleton test:", id(piCnst), id(piCnst1), id(piCnst)-id(piCnst1))
if piCnst != piCnst1:
    raise Exception("PionConst is not a singleton class!")

##! Check built-in methods:
PionConstTest = 2
print()
print("PionConstTest:", PionConstTest, " check built-in methods.")


#  redirect standard out to a file
restoreOut = sys.stdout
sys.stdout = open("Scratch/pionTst.out","w")
print(piCnst)
sys.stdout.close()
sys.stdout = restoreOut
# compare the standard output to a reference file
a = os.popen('diff Scratch/pionTst.out 02-Tests/pionTst.ref')
output = a.read()
if (output == ""):
    pass
else:
    testFails = testFails + 1

##! Check get methods:
PionConstTest = 3
print()
print("PionConstTest:", PionConstTest, " check get methods.")
print("----> print() method; tests all get methods")
piCnst.print()

##! Complete:
print()
print("========  PionConst: tests complete  ========")
if (testFails == 0):
  sys.exit(0)
else:
  sys.exit(1)
