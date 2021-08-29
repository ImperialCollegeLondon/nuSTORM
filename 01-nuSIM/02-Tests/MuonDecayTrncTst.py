#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "MuonDecay" class -- for truncated decays
=================================

  Assumes that nuSim code is in python path.

  Script starts by testing built in methods.  Then a soak test with a
  large number of decays is executed.  Finally a set of reference plots
  are generated.

"""

import MuonConst as muCnst
import MuonDecay as md
import Simulation as Simu
import matplotlib.pyplot as plt
import numpy as np
import math as mt

##! Start:
print("========  MuonDecayTruncated: tests start  ========")

mc = muCnst.MuonConst()

##! Create instance, test built-in methods:
MuonDecayTruncatedTest = 1
print()
print("MuonDecayTruncatedTest:", MuonDecayTruncatedTest, " Create muon decay, print quantities.")
print("    ----> no Tmax argument")
Dcy=md.MuonDecay()
print("    __str__:", Dcy)
print("    --repr__", repr(Dcy))
del Dcy

Tmax = 2.*mc.lifetime()
print("    ----> Tmax =", Tmax, "ns")
Dcy=md.MuonDecay(Tmax=Tmax)
print("    __str__:", Dcy)
print("    --repr__", repr(Dcy))
del Dcy

##! Create instance, test dynamic methods:
MuonDecayTruncatedTest = 2
print()
print("MuonDecayTruncatedTest:", MuonDecayTruncatedTest, " Create muon decay, print quantities.")
Dcy=md.MuonDecay(Tmax=Tmax)
print("     _Lifetime:", Dcy.getLifetime())
print("          _v_e:", Dcy._v_e)
print("        _v_nue:", Dcy._v_nue)
print("       _v_numu:", Dcy._v_numu)
SumE = Dcy._v_e[0]+Dcy._v_nue[0]+Dcy._v_numu[0]
DifE = SumE - mc.mass()
print("    Sum energy:", SumE, "; difference to muon mass:", DifE)
del Dcy

##! Soak test, generate many decays:
MuonDecayTruncatedTest = 3
print()
print("MuonDecayTruncatedTest:", MuonDecayTruncatedTest, " Create many decays.")
Dcy = []
for i in range(100000):
    Dcy.append(md.MuonDecay(Tmax=Tmax))
for i in range(5):
    print(Dcy[i])

##! Plot result of soak test:
MuonDecayTruncatedTest = 4
print()
print("MuonDecayTruncatedTest:", MuonDecayTruncatedTest, " plots from soak test.")

t       = np.array([])
Ee      = np.array([])
Enue    = np.array([])
Enumu   = np.array([])
cosnue  = np.array([])
cosnumu = np.array([])
s = 0.
for muDcy in Dcy:
    t     = np.append(t,     muDcy.getLifetime())
    Ee    = np.append(Ee,    muDcy.get4ve()[0])
    Enue  = np.append(Enue,  muDcy.get4vnue()[0])
    Enumu = np.append(Enumu, muDcy.get4vnumu()[0])
    
    p_e    = muDcy.get4ve()[1]
    p_nue  = muDcy.get4vnue()[1]
    p_numu = muDcy.get4vnumu()[1]

    mag_e    = np.linalg.norm(p_e)
    mag_nue  = np.linalg.norm(p_nue)
    mag_numu = np.linalg.norm(p_numu)

    p_e    = p_e    / mag_e
    p_nue  = p_nue  / mag_nue
    p_numu = p_numu / mag_numu
                    
    cos_e   = np.dot(p_e, p_nue)
    cos_mu  = np.dot(p_e, p_numu)
    cosnue  = np.append(cosnue,  cos_e)
    cosnumu = np.append(cosnumu, cos_mu)
    s += 1.

#-- Lifetime distribution:
n, bins, patches = plt.hist(t, bins=50, color='y', range=(0.,2.*Tmax), log=True)
plt.xlabel('Time (s)')
plt.ylabel('Entries')
plt.title('Lifetime distribution')
# add a 'best fit' line
l = 1./mc.lifetime()
y = n[0]*np.exp(-l*bins)
plt.plot(bins, y, '-', color='b')
plt.savefig('Scratch/MuonDecayTrncTst_plot1.pdf')
plt.close()

#-- Energy distributions:
n, bins, patches = plt.hist(Ee, bins=50, color='y', range=(0.,50.0))
plt.xlabel('Energy (MeV)')
plt.ylabel('Frequency')
plt.title('Electron energy distribution')
# add a 'best fit' line
x2 = bins*bins
x3 = x2*bins
y  = 3.*x2 - (4./mc.mass())*x3
y  = (16.*s/(mc.mass()**3)) * y
plt.plot(bins, y, '-', color='b')
plt.savefig('Scratch/MuonDecayTrncTst_plot2.pdf')
plt.close()

n, bins, patches = plt.hist(Enue, bins=50, color='y', range=(0.,50.))
plt.xlabel('Energy (MeV)')
plt.ylabel('Frequency')
plt.title('Electron-neutrino energy distribution')
# add a 'best fit' line
y  = x2 - (2./mc.mass())*x3
y  = (96.*s/(mc.mass()**3)) * y
plt.plot(bins, y, '-', color='b')
plt.savefig('Scratch/MuonDecayTrncTst_plot3.pdf')
plt.close()

n, bins, patches = plt.hist(Enumu, bins=50, color='y', range=(0.,50.))
plt.xlabel('Energy (MeV)')
plt.ylabel('Frequency')
plt.title('Muon-neutrino energy distribution')
# add a 'best fit' line
y  = 3.*x2 - (4./mc.mass())*x3
y  = (16.*s/(mc.mass()**3)) * y
plt.plot(bins, y, '-', color='b')
plt.savefig('Scratch/MuonDecayTrncTst_plot4.pdf')
plt.close()

#-- Angular distributions:
n, bins, patches = plt.hist(cosnue, bins=50, color='y', range=(-1.,1.))
plt.xlabel('Cos(theta)')
plt.ylabel('Frequency')
plt.title('Electron neutrino cos(theta) distribution')
plt.savefig('Scratch/MuonDecayTrncTst_plot5.pdf')
plt.close()

n, bins, patches = plt.hist(cosnumu, bins=50, color='y', range=(-1.,1.))
plt.xlabel('Cos(theta)')
plt.ylabel('Frequency')
plt.title('Muon neutrino cos(theta) distribution')
plt.savefig('Scratch/MuonDecayTrncTst_plot6.pdf')
plt.close()

##! Complete:
print()
print("========  MuonDecayTruncated: tests complete  ========")
