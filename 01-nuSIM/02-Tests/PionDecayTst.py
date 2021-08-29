#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "PionDecay" class
=================================

  Assumes that nuSim code is in python path.

  Script starts by testing built in methods.  Then a soak test with a
  large number of decays is executed.  Finally a set of reference plots
  are generated.

"""

import PionConst as piCnst
import PionDecay as pd
import Simulation as Simu
import matplotlib.pyplot as plt
import numpy as np
import math as mt

##! Start:
print("========  PionDecay: tests start  ========")

##! Create instance, test built-in methods:
PionDecayTest = 1
print()
print("PionDecayTest:", PionDecayTest, " Create pion decay, print quantities.")
Dcy=pd.PionDecay()
print("    __str__:", Dcy)
print("    --repr__", repr(Dcy))
del Dcy

##! Create instance, test dynamic methods:
PionDecayTest = 2
print()
print("PionDecayTest:", PionDecayTest, " Create pion decay, print quantities.")
Dcy=pd.PionDecay()
print("     _Lifetime:", Dcy.getLifetime())
print("         _v_mu:", Dcy._v_mu)
print("       _v_numu:", Dcy._v_numu)
pc = piCnst.PionConst()
SumE = Dcy._v_mu[0] + Dcy._v_numu[0]
DifE = SumE - pc.mass()
print("    Sum energy:", SumE, "; difference to muon mass:", DifE)
del Dcy

##! Soak test, generate many decays:
PionDecayTest = 3
print()
print("PionDecayTest:", PionDecayTest, " Create many decays.")
Dcy = []
for i in range(50000):
    Dcy.append(pd.PionDecay())
for i in range(5):
    print(Dcy[i])

##! Plot result of soak test:
PionDecayTest = 4
print()
print("PionDecayTest:", PionDecayTest, " plots from soak test.")

t         = np.array([])
Emu       = np.array([])
Enumu     = np.array([])
cosTheta  = np.array([])
phi       = np.array([])
s = 0.
for piDcy in Dcy:
    t      = np.append(t,     piDcy.getLifetime())
    Emu    = np.append(Emu,   piDcy.get4vmu()[0])
    Enumu  = np.append(Enumu, piDcy.get4vnumu()[0])

    p_mu   = piDcy.get4vmu()[1]
    p_numu = piDcy.get4vnumu()[1]

    mag_mu    = np.linalg.norm(p_mu)
    
    phi = np.append(phi,     mt.atan2(p_mu[1],p_mu[0]))
    cosTheta = np.append(cosTheta,     p_mu[2]/mag_mu)
                    
#-- Lifetime distribution:
n, bins, patches = plt.hist(t, bins=50, color='y', range=(0.,2.E-8), log=True)
plt.xlabel('Time (s)')
plt.ylabel('Entries')
plt.title('Lifetime distribution')
# add a 'best fit' line
l = 1./pc.lifetime()
y = n[0]*np.exp(-l*bins)
plt.plot(bins, y, '-', color='b')
plt.show()

n, bins, patches = plt.hist(Emu, bins=50, color='y', range=(100.,120.))
plt.xlabel('Energy (MeV)')
plt.ylabel('Frequency')
plt.title('muon energy distribution')
plt.show()

n, bins, patches = plt.hist(Enumu, bins=50, color='y', range=(20.,40.))
plt.xlabel('Energy (MeV)')
plt.ylabel('Frequency')
plt.title('Muon-neutrino energy distribution')
plt.show()

#-- Angular distributions:
n, bins, patches = plt.hist(phi, bins=50, color='y', range=(-4.,4.))
plt.xlabel('Phi')
plt.ylabel('Frequency')
plt.title('muon phi distribution')
plt.show()

n, bins, patches = plt.hist(cosTheta, bins=50, color='y', range=(-1.1,1.1))
plt.xlabel('Cos(theta)')
plt.ylabel('Frequency')
plt.title('Muon cos(theta) distribution')
plt.show()

##! Complete:
print()
print("========  PionDecay: tests complete  ========")
