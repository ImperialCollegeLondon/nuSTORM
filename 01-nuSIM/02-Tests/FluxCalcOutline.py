#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "FluxShape"
===========================

  Assumes that nuSim code is in path.

"""

import os
import MuonConst as muCnst
import NeutrinoEventInstance as nuEvtInst
import numpy as np
import matplotlib.pyplot as plt
import Simulation as Simu

mc = muCnst.MuonConst()

##! Start:
print("========  FluxShape: tests start  ========")

##! Create instance, test built-in methods:
FluxShapeTest = 1
print()
print("FluxShapeTest:", FluxShapeTest, " get one event.")
NoEvtAccepted = True
Prt1True = False
Prt2True = False

nuSIMPATH = os.getenv('nuSIMPATH')
filename  = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')

# .. assume mu+
muBmE   = 3.8
nueSlp  = 1./(muBmE - 0.)
numuSlp  = 1./(muBmE - mc.mass()/1000.)

while NoEvtAccepted:
    nuEI    = nuEvtInst.NeutrinoEventInstance(muBmE, filename)
    pt_nue  = np.sqrt(nuEI.getnue4mmtm()[1][0]**2 + nuEI.getnue4mmtm()[1][1]**2)
    pt_numu = np.sqrt(nuEI.getnumu4mmtm()[1][0]**2 + nuEI.getnumu4mmtm()[1][1]**2)

    Rr      = (250.-nuEI.getTraceSpaceCoord()[0]+50.)*(pt_nue/nuEI.getnue4mmtm()[1][2])
    if abs(Rr) < 5.:
        Rn1  = Simu.getRandom()
        Tst1 = nueSlp*nuEI.getnue4mmtm()[0]
        if Rn1 < Tst1:
            print("    Electron! Rn1, Tst1, nuEI:", Rn1, Tst1, "\n", nuEI.__str__())
            Prt1True = True

    Rr      = (250.-nuEI.getTraceSpaceCoord()[0]+50.)*(pt_numu/nuEI.getnumu4mmtm()[1][2])
    if abs(Rr) < 5.:      
        Rn1  = Simu.getRandom()
        Tst2 = numuSlp*(nuEI.getnumu4mmtm()[0] - mc.mass()/1000.)
        if Rn1 < Tst2:
            print("    Muon! Rn1, Tst2, nuEI:", Rn1, Tst2, "\n", nuEI.__str__())
            Prt2True = True
            if Prt2True and Prt1True:
                NoEvtAccepted = False
    del nuEI
        
##! Test 
FluxShapeTest = 2
print()
print("FluxShapeTest:", FluxShapeTest, \
      "soak test.")
nuEIelec = []
nuEImuon = []
for i in range(10000):
    nuEI = nuEvtInst.NeutrinoEventInstance(muBmE, filename)
    pt_nue  = np.sqrt(nuEI.getnue4mmtm()[1][0]**2 + nuEI.getnue4mmtm()[1][1]**2)
    pt_numu = np.sqrt(nuEI.getnumu4mmtm()[1][0]**2 + nuEI.getnumu4mmtm()[1][1]**2)

    Rr      = (250.-nuEI.getTraceSpaceCoord()[0]+50.)*(pt_nue/nuEI.getnue4mmtm()[1][2])
    if abs(Rr) < 5.:
        Rn1  = Simu.getRandom()
        Tst1 = nueSlp*nuEI.getnue4mmtm()[0]
        if Rn1 < Tst1:
            nuEIelec = np.append(nuEIelec, nuEI.getnue4mmtm()[0])

    Rr      = (250.-nuEI.getTraceSpaceCoord()[0]+50.)*(pt_numu/nuEI.getnumu4mmtm()[1][2])
    if abs(Rr) < 5.:      
        Rn1  = Simu.getRandom()
        Tst2 = numuSlp*(nuEI.getnumu4mmtm()[0] - mc.mass()/1000.)
        if Rn1 < Tst2:
            nuEImuon = np.append(nuEImuon, nuEI.getnumu4mmtm()[0])

#-- Energy distributions:
n, bins, patches = plt.hist(nuEIelec, bins=50, color='y', range=(0.,5.))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('Electron-neutrino energy distribution')
plt.savefig('Scratch/FluxCalcOutline_plot1.pdf')
plt.close()

n, bins, patches = plt.hist(nuEImuon, bins=50, color='y', range=(0.,5.))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('Muon-neutrino energy distribution')
plt.savefig('Scratch/FluxCalcOutline_plot2.pdf')
plt.close()

##! Complete:
print()
print("========  FluxShape: tests complete  ========")
