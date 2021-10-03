#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "Simulation" class ... initialisation and get methods
==================================

  Simulation.py -- set "relative" path to code

"""

import os
import Simulation as Simu
import nuSTORMPrdStrght as nuStrt

##! Start:
print("========  Simulation: tests start  ========")

##! Test singleton class feature:
SimulationTest = 1
print()
print("SimulationTest:", SimulationTest, " check if class is a singleton.")
nuSIMPATH    = os.getenv('nuSIMPATH')
filename     = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
rootfilename = os.path.join(nuSIMPATH, 'Scratch/nuSIM-Simulation-tst.root')
print(filename)
print(rootfilename)
Smltn  = Simu.Simulation(100, 5., filename, rootfilename)
Smltn1 = Simu.Simulation(100, 5., filename, rootfilename)
print("---->Smltn singleton test:", id(Smltn), id(Smltn1), id(Smltn)-id(Smltn1))
if Smltn != Smltn1:
    raise Exception("Simulation is not a singleton class!")

##! Check get methods:
SimulationTest = 2
print()
print("SimulationTest:", SimulationTest, " check get methods.")
print("    Simulation: version:", Smltn.CdVrsn())
rS = Smltn.getRandomSeed()
if (rS < 100000000.0):
      raise Exception("Simulation seed is too small ", rS)
print("    __RandomSeed:", rS)
r1 = Simu.getRandom()
if ((r1 < 0.0) or (r1 > 1.0)):
    raise Exception("Simulation returns random number not beteen 0.0 and 1.0 ", r1)
print("    Random no.  :", r1)
r1 = Simu.getRandom()
if ((r1 < 0.0) or (r1 > 1.0)):
    raise Exception("Simulation returns random number not beteen 0.0 and 1.0 ", r1)
print("    Random no.  :", r1)
r1 = Simu.getRandom()
if ((r1 < 0.0) or (r1 > 1.0)):
    raise Exception("Simulation returns random number not beteen 0.0 and 1.0 ", r1)
print("    Random no.  :", r1)
r1 = Simu.getRandom()
if ((r1 < 0.0) or (r1 > 1.0)):
    raise Exception("Simulation returns random number not beteen 0.0 and 1.0 ", r1)
print("    Random no.  :", r1)
r1 = Simu.getRandom()
if ((r1 < 0.0) or (r1 > 1.0)):
    raise Exception("Simulation returns random number not beteen 0.0 and 1.0 ", r1)
print("    Random no.  :", r1)
r1 = Simu.getRandom()
if ((r1 < 0.0) or (r1 > 1.0)):
    raise Exception("Simulation returns random number not beteen 0.0 and 1.0 ", r1)
print("    Random no.  :", r1)
r1 = Simu.getRandom()
if ((r1 < 0.0) or (r1 > 1.0)):
    raise Exception("Simulation returns random number not beteen 0.0 and 1.0 ", r1)
print("    Random no.  :", r1)
r1 = Simu.getRandom()
if ((r1 < 0.0) or (r1 > 1.0)):
    raise Exception("Simulation returns random number not beteen 0.0 and 1.0 ", r1)
print("    Random no.  :", r1)
r1 = Simu.getRandom()
if ((r1 < 0.0) or (r1 > 1.0)):
    raise Exception("Simulation returns random number not beteen 0.0 and 1.0 ", r1)
print("    Random no.  :", r1)
r1 = Simu.getRandom()
if ((r1 < 0.0) or (r1 > 1.0)):
    raise Exception("Simulation returns random number not beteen 0.0 and 1.0 ", r1)
print("    Random no.  :", r1)

##! Complete:
print()
print("========  Simulation: tests complete  ========")
