#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "nuSTORMConst" class ... initialisation and get methods
=================================

  nuSTORMConst.py -- set "relative" path to code

"""

import nuSTORMConst as nc
import sys, os

testFails = 0
##! Start:
print("========  nuSTORMConst: tests start  ========")


##! Test singleton class feature:
nuSTORMConstTest = 1
print()
print("nuSTORMConstTest:", nuSTORMConstTest, " check if class is a singleton.")
nuSTRMCnst  = nc.nuSTORMConst()
nuSTRMCnst1 = nc.nuSTORMConst()
print("    nuSTRMCnst singleton test:", id(nuSTRMCnst), id(nuSTRMCnst1), id(nuSTRMCnst)-id(nuSTRMCnst1))
if nuSTRMCnst != nuSTRMCnst1:
    raise Exception("nuSTORMConst is not a singleton class!")

##! Check built-in methods:
nuSTORMConstTest = 2
print()
print("nuSTORMConstTest:", nuSTORMConstTest, " check built-in methods.")


#  redirect standard out to a file
#restoreOut = sys.stdout
#sys.stdout = open("Scratch/pionTst.out","w")
#print(nuSTRMCnst)
#sys.stdout.close()
#sys.stdout = restoreOut
## compare the standard output to a reference file
#a = os.popen('diff Scratch/pionTst.out 02-Tests/referenceOutput/pionTst.ref')
#output = a.read()
#if (output == ""):
#    pass
#else:
#    testFails = testFails + 1

##! Check get methods:
nuSTORMConstTest = 3
print()
print("nuSTORMConstTest:", nuSTORMConstTest, " check get methods.")
print("----> print() method; tests all get methods")
print("")
nuSTRMCnst.print()

##! Check get pdgCode:
#nuSTORMConstTest = 5
#if (nuSTRMCnst.pdgCode() != 211):
#  testFails = testFails + 1


##! Complete:
print()
print("========  nuSTORMConst: tests complete  ========")
print ("testFails is ", testFails)
if (testFails == 0):
  sys.exit(0)
else:
  sys.exit(1)
