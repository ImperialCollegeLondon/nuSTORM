#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "MuonBeam4CoolingDemo" class
============================================

  Assumes python path includes nuSim code.

  Script starts by testing built in methods.  Then a soak test with a set 
  of reference plots.

"""

import numpy as np
import matplotlib.pyplot as plt
import MuonBeam4CoolingDemo as muBm4ClDmo
import os

##! Start:
print("========  MuonBeam4CoolingDemo: tests start  ========")

##! Test singleton class feature:
MuonBeam4CoolingDemoTest = 1
print()
print("MuonBeam4CoolingDemoTest:", MuonBeam4CoolingDemoTest, " check if class is a singleton.")
nuSIMPATH = os.getenv('nuSIMPATH')
print(nuSIMPATH)
filename  = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-Muons4CoolingDemo-Params-v1.0.csv')
mu4CD  = muBm4ClDmo.MuonBeam4CoolingDemo(filename)
mu4CD1 = muBm4ClDmo.MuonBeam4CoolingDemo(filename)
print("    mu4CD singleton test: result should be 0:", id(mu4CD)-id(mu4CD1))
if mu4CD != mu4CD1:
    raise Exception("MuonBeam4CoolingDemo is not a singleton class!")

##! Check built-in methods:
MuonBeam4CoolingDemoTest = 2
print()
print("MuonBeam4CoolingDemoTest:", MuonBeam4CoolingDemoTest, " check built-in methods.")
print("    __repr__:", repr(mu4CD))
print("    __str__:")
print(mu4CD)

##! Test muon instance generator
MuonBeam4CoolingDemoTest = 3
import mu4ClDmo as mu4ClDmo
print()
print("MuonBeam4CoolingDemoTest:", MuonBeam4CoolingDemoTest, " test instanciation of muon for rest beam.")
muInst = mu4ClDmo.mu4ClDmo()

##! Check built-in methods:
MuonBeam4CoolingDemoTest = 4
print()
print("MuonBeam4CoolingDemoTest:", MuonBeam4CoolingDemoTest, " check built-in methods of muon instance generator.")
print("    __repr__:", repr(muInst))
print("    __str__:")
print(muInst)
del muInst

##! Soak test:
MuonBeam4CoolingDemoTest = 5
print()
print("MuonBeam4CoolingDemoTest:", MuonBeam4CoolingDemoTest, " soak test; generate muon beam.")
muBm = []
for i in range(100000):
    muBm.append(mu4ClDmo.mu4ClDmo())
for i in range(5):
    print(muBm[i])

##! Plot result of soak test:
MuonBeam4CoolingDemoTest = 5
print()
print("MuonBeam4CoolingDemoTest:", MuonBeam4CoolingDemoTest, " plot result of soak test.")

tGen    = np.array([])
xGen    = np.array([])
yGen    = np.array([])
pxGen   = np.array([])
pyGen   = np.array([])
pzGen   = np.array([])
EGen    = np.array([])

s = 0.
for muon in muBm:
    Rmu   = muon.Rmu()
    Pmu   = muon.Pmu()
    
    tGen  = np.append(tGen,   Rmu[0])
    xGen  = np.append(xGen,   Rmu[1][0])
    yGen  = np.append(yGen,   Rmu[1][1])

    EGen  = np.append(EGen,   Pmu[0])
    pxGen = np.append(pxGen,  Pmu[1][0])
    pyGen = np.append(pyGen,  Pmu[1][1])
    pzGen = np.append(pzGen,  Pmu[1][2])

    s += 1.

#-- Lifetime distribution:
n, bins, patches = plt.hist(tGen, bins=50, color='y', range=(-15.,15.), log=False)
plt.xlabel('t (ns)')
plt.ylabel('Entries')
plt.title('Time distribution')
plt.savefig('Scratch/MuonBeam4CoolingDemoTst_plot1.pdf')
plt.close()

n, bins, patches = plt.hist(xGen, bins=50, color='y', range=(-20.,20.), log=False)
plt.xlabel('x (mm)')
plt.ylabel('Entries')
plt.title('x position distribution')
plt.savefig('Scratch/MuonBeam4CoolingDemoTst_plot2.pdf')
plt.close()

n, bins, patches = plt.hist(yGen, bins=50, color='y', range=(-20.,20.), log=False)
plt.xlabel('y (mm)')
plt.ylabel('Entries')
plt.title('y position distribution')
plt.savefig('Scratch/MuonBeam4CoolingDemoTst_plot3.pdf')
plt.close()

n, bins, patches = plt.hist(EGen, bins=50, color='y', range=(120,170.), log=False)
plt.xlabel('E (MeV)')
plt.ylabel('Entries')
plt.title('Energy distribution')
plt.savefig('Scratch/MuonBeam4CoolingDemoTst_plot4.pdf')
plt.close()

n, bins, patches = plt.hist(pxGen, bins=50, color='y', range=(-30.,30.), log=False)
plt.xlabel('px (MeV/c)')
plt.ylabel('Entries')
plt.title('px distribution')
plt.savefig('Scratch/MuonBeam4CoolingDemoTst_plot5.pdf')
plt.close()

n, bins, patches = plt.hist(pyGen, bins=50, color='y', range=(-30.,30.), log=False)
plt.xlabel('py (MeV/c)')
plt.ylabel('Entries')
plt.title('py distribution')
plt.savefig('Scratch/MuonBeam4CoolingDemoTst_plot6.pdf')
plt.close()

n, bins, patches = plt.hist(pzGen, bins=50, color='y', range=(70,130.), log=False)
plt.xlabel('pz (MeV/c)')
plt.ylabel('Entries')
plt.title('pz distribution')
plt.savefig('Scratch/MuonBeam4CoolingDemoTst_plot7.pdf')
plt.close()



##! Complete:
print()
print("========  MuonBeam4CoolingDemo: tests complete  ========")


