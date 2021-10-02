#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "PionConst" class ... initialisation and get methods
=================================

  PionConst.py -- set "relative" path to code

"""

import PionConst as pc

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
print("    __repr__:")
print(piCnst)

##! Check get methods:
PionConstTest = 3
print()
print("PionConstTest:", PionConstTest, " check get methods.")
print("----> print() method; tests all get methods")
piCnst.print()

##! Complete:
print()
print("========  PionConst: tests complete  ========")
