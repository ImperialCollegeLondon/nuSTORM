#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "PionEventInstance" class
=============================================

  Assumes that nuSim code is in path.

  Script starts by testing built in methods.  Then a soak test with a
  large number of decays is executed.  Finally a set of reference plots
  are generated.

"""

import MuonConst as muCnst
import PionConst as piCnst
import PionEventInstance as piEvtInst
import numpy as np
import matplotlib.pyplot as plt
import Simulation as Simu

mc = muCnst.MuonConst()
pc = piCnst.PionConst()

##! Start:
print("========  PionEventInstance: tests start  ========")

##! Create instance, test built-in methods:
PionEventInstanceTest = 1
print()
print("PionEventInstanceTest:", PionEventInstanceTest, " Test built-in methods.")
piEI = piEvtInst.PionEventInstance(8.)
print("    __str__:", piEI)
print("    --repr__", repr(piEI))

##! Test get/set methods:
PionEventInstanceTest = PionEventInstanceTest + 1
print()
print("PionEventInstanceTest:", PionEventInstanceTest, \
      "Test get/set methods.")
print("    Pion momentum (GeV):", piEI.getppi())

##! Test methods by which neutrino-creation event is generated:
PionEventInstanceTest = PionEventInstanceTest + 1
print()
print("PionEventInstanceTest:", PionEventInstanceTest, \
      "Test methods by which muon-creation event is generated.")
decayPnt, PPi, Pmu, Pnumu, tpi = piEI.CreateMuon()
print("    Neutrino event: trace-space coordinates of muon at decay, P_pi, P_mu, P_numu:", decayPnt, "\n", \
           "    ", PPi,"    ", Pmu, "    ", Pnumu)
ETot = Pmu[0]+Pnumu[0]
PTot = Pmu[1]+Pnumu[1]
print ("totol energy is ", ETot)
print ("total momentum is ", PTot)
del piEI

##! Soak test:
PionEventInstanceTest = PionEventInstanceTest + 1
print()
print("PionEventInstanceTest:", PionEventInstanceTest, \
      "soak test.")
piEI = []
Ppion=6.0
for i in range(10000):
    piEI.append(piEvtInst.PionEventInstance(Ppion))
for i in range(5):
    print("    piEI[i]:", piEI[i])

#! Plot result of soak test:
PionEventInstanceTest = PionEventInstanceTest + 1
NeutrinoEventInstanceTest = 5
print()
print("PionEventInstanceTest:", PionEventInstanceTest, \
      " plots from soak test.")

tpi      = np.array([])
s        = np.array([])
Epi      = np.array([])
Emu      = np.array([])
Enumu    = np.array([])
phi      = np.array([])
costheta = np.array([])


for piEvt in piEI:
    tpi   = np.append(tpi,    piEvt.gettpi())
    s     = np.append(s,     piEvt.getTraceSpaceCoord()[0])
#    Ee    = np.append(Ee,    nuEvt.gete4mmtm()[0])
#    Enue  = np.append(Enue,  nuEvt.getnue4mmtm()[0])
    Emu = np.append(Emu, piEvt.getmu4mmtm()[0])
    Enumu = np.append(Enumu, piEvt.getnumu4mmtm()[0])
    phi = np.append(phi, piEvt.getphi())
    costheta = np.append(costheta, piEvt.getcostheta())
    
#    pt_e    = np.sqrt(nuEvt.gete4mmtm()[1][0]**2 + nuEvt.gete4mmtm()[1][1]**2)
#    pt_nue  = np.sqrt(nuEvt.getnue4mmtm()[1][0]**2 + nuEvt.getnue4mmtm()[1][1]**2)
#    pt_numu = np.sqrt(nuEvt.getnumu4mmtm()[1][0]**2 + nuEvt.getnumu4mmtm()[1][1]**2)

#    tane    = np.append(tane,    (pt_e/nuEvt.gete4mmtm()[1][2]) )
#    tannue  = np.append(tannue,  (pt_nue/nuEvt.getnue4mmtm()[1][2]) )
#    tannumu = np.append(tannumu, (pt_numu/nuEvt.getnumu4mmtm()[1][2]) )

#-- Lifetime distribution:
n, bins, patches = plt.hist(s, bins=50, color='y', log=True)
plt.xlabel('s (m)')
plt.ylabel('Entries')
plt.title('s distribution')


## add a 'best fit' line - convert en
charDist = Ppion*1000.0*pc.SoL()*pc.lifetime()/pc.mass()
slope = 1./charDist
y = n[0]*np.exp(-slope*bins)
plt.plot(bins, y, '-', color='b')
plt.show()
plt.close()

n, bins, patches = plt.hist(tpi, bins=50, color='y', log=True)
plt.xlabel('tpi (ns)')
plt.ylabel('Entries')
plt.title('Time of pion decay distribution')

## add a 'best fit' line - convert en
charDist = pc.lifetime()
slope = 1./charDist
y = n[0]*np.exp(-slope*bins)
plt.plot(bins, y, '-', color='b')
plt.show()
plt.close()

#-- Energy distributions:
n, bins, patches = plt.hist(Emu, bins=50, color='y', range=(0.,10.0))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('muon energy distribution')
plt.show()
plt.close()

n, bins, patches = plt.hist(Enumu, bins=50, color='y', range=(0.,10.))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('Muon-neutrino energy distribution')
plt.show()
plt.close()

n, bins, patches = plt.hist(phi, bins=50, color='y', range=(0.,6.3))
plt.xlabel('Angle')
plt.ylabel('Frequency')
plt.title('Muon phi distribution')
plt.show()
plt.close()

n, bins, patches = plt.hist(costheta, bins=50, color='y', range=(-2.,2.))
plt.xlabel('Angle')
plt.ylabel('Frequency')
plt.title('Muon phi distribution')
plt.show()
plt.close()

plt.scatter(Emu,phi)
plt.xlabel('Emu (GeV)')
plt.ylabel('phi')
plt.title('phi vs Emu plot')
plt.show()
plt.close()
##! Complete:
print()
print("========  PionEventInstance: tests complete  ========")
