#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "MuonConst" class ... initialisation and get methods
=================================

  MuonConst.py -- set "relative" path to code

"""

import MuonConst as mc
import sys

##! Start:
print("========  MuonConst: tests start  ========")

##! Test singleton class feature:
MuonConstTest = 1
print()
print("MuonConstTest:", MuonConstTest, " check if class is a singleton.")
muCnst  = mc.MuonConst()
muCnst1 = mc.MuonConst()
print("    muCnst singleton test:", id(muCnst), id(muCnst1), id(muCnst)-id(muCnst1))
if muCnst != muCnst1:
    raise Exception("MuonConst is not a singleton class!")

##! Check built-in methods:
MuonConstTest = 2
print()
print("MuonConstTest:", MuonConstTest, " check built-in methods.")
print("    __repr__:")
print(muCnst)

##! Check get methods:
MuonConstTest = 3
print()
print("MuonConstTest:", MuonConstTest, " check get methods.")
print("----> print() method; tests all get methods")
muCnst.print()

##! Check get pdgCode:
MuonConstTest = 4
if (muCnst.pdgCode() != 13):
  sys.exit("wrong pdgcode")

##! Complete:
print()
print("========  MuonConst: tests complete  ========")
