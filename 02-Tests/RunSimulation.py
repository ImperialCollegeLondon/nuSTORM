#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "Simulation" class ... initialisation and get methods
==================================

  Simulation.py -- set "relative" path to code

"""

import os
import Simulation as Simu

##! Start:
print("========  Simulation: start  ========")
print()

nuSIMPATH = os.getenv('nuSIMPATH')
filename  = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
rootfilename = os.path.join(nuSIMPATH, 'Scratch/nuSIM-RunSimulation.root')

Smltn = Simu.Simulation(5000, 6., filename, rootfilename)

print()
print(".... Execute simulation")

Smltn.RunSim()

##! Complete:
print()
print("========  Simulation: complete  ========")
