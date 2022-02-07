#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "PionEventInstance" class
=============================================

  Assumes that nuSim code is in path.

  Script starts by testing built in methods.  Then a soak test with a
  large number of decays is executed.  Finally a set of reference plots
  are generated.

  Version 1.1            Paul Kyberd                    20 October 2021
  Version 1.2            Marvin Pfaff                   20 January 2022

  Add sys.exit("message ") to test 4

"""

import MuonConst as muCnst
import PionConst as piCnst
import PionEventInstance as piEvtInst
import nuSTORMPrdStrght as nuPrdStrt
import nuSTORMTrfLineCmplx as nuTrfLineCmplx
import particle as particle
import eventHistory
import RandomGenerator as Rndm
import numpy as np
import matplotlib.pyplot as plt
import Simulation as Simu
import sys, os

mc = muCnst.MuonConst()
pc = piCnst.PionConst()
nuStrt = nuPrdStrt.nuSTORMPrdStrght('11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
nuSIMPATH = os.getenv('nuSIMPATH')
trfCmplxFile = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-TrfLineCmplx-Params-v1.0.csv')
nuTrLnCmplx = nuTrfLineCmplx.nuSTORMTrfLineCmplx(trfCmplxFile)

##! Start:
print("========  PionEventInstance: tests start  ========")

n_events = 5000

##! Create instance, test built-in methods:
PionEventInstanceTest = 1
print()
print("PionEventInstanceTest:", PionEventInstanceTest, " Test built-in methods.")
piMomentum = 5.
piEI = piEvtInst.PionEventInstance(piMomentum)
print("    __str__:", piEI)
print("    --repr__", repr(piEI))

##! Test get/set methods:
PionEventInstanceTest = PionEventInstanceTest + 1
print()
print("PionEventInstanceTest:", PionEventInstanceTest, \
      "Test get/set methods.")
print("    Pion momentum (GeV):", piEI.getppi())
if (piEI.getppi() != piMomentum):
  sys.exit("getppi returns the wrong value")


##! Test methods by which muon-creation event is generated:
PionEventInstanceTest = PionEventInstanceTest + 1
print()
print("PionEventInstanceTest:", PionEventInstanceTest, \
      "Test methods by which muon-creation event is generated.")
decayPnt, PPi, Pmu, Pnumu = piEI.CreateMuon()
print("    Muon event: trace-space coordinates of pion at decay, P_pi, P_mu, P_numu:", decayPnt, "\n", \
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
Ppion=piMomentum
for i in range(n_events):
    piEI.append(piEvtInst.PionEventInstance(Ppion))
for i in range(5):
    print("    piEI[i]:", piEI[i])

#! Plot result of soak test:
PionEventInstanceTest = PionEventInstanceTest + 1
print()
print("PionEventInstanceTest:", PionEventInstanceTest, \
      " plots from soak test.")

s        = np.array([])
pxpi      = np.array([])
pypi      = np.array([])
pzpi      = np.array([])
Emu      = np.array([])
pxmu      = np.array([])
pymu      = np.array([])
pzmu      = np.array([])
Enumu    = np.array([])
pxnumu      = np.array([])
pynumu      = np.array([])
pznumu      = np.array([])
phi      = np.array([])
costheta = np.array([])


for piEvt in piEI:
    s     = np.append(s,     piEvt.getTraceSpaceCoord()[0])
    pxpi = np.append(pxpi,   piEvt.getTraceSpaceCoord()[4]*piEvt.getppi())
    pypi = np.append(pypi,   piEvt.getTraceSpaceCoord()[5]*piEvt.getppi())
    pzpi = np.append(pzpi,   np.sqrt(piEvt.getppiGen()**2 - pxpi**2 - pypi**2))
    Emu = np.append(Emu, piEvt.getmu4mmtm()[0])
    pxmu = np.append(pxmu, piEvt.getmu4mmtm()[1][0])
    pymu = np.append(pymu, piEvt.getmu4mmtm()[1][1])
    pzmu = np.append(pzmu, piEvt.getmu4mmtm()[1][2])
    Enumu = np.append(Enumu, piEvt.getnumu4mmtm()[0])
    pxnumu = np.append(pxnumu, piEvt.getnumu4mmtm()[1][0])
    pynumu = np.append(pynumu, piEvt.getnumu4mmtm()[1][1])
    pznumu = np.append(pznumu, piEvt.getnumu4mmtm()[1][2])
    phi = np.append(phi, piEvt.getphi())
    costheta = np.append(costheta, piEvt.getcostheta())


##! Soak test with distributions from histogram:
PionEventInstanceTest = PionEventInstanceTest + 1
print()
print("PionEventInstanceTest:", PionEventInstanceTest, \
      "soak test generating distributions from particle.")
piEI_h = []
Ppion=piMomentum
##Use the following 3 lines to create according to histogram distributions
#nuSIMPATH = os.getenv('nuSIMPATH')
#rootInputFilename = os.path.join(nuSIMPATH, 'Scratch/plots_endProdStrght.root')
#RndmGen = Rndm.RandomGenerator(rootInputFilename,'histP','histXPS')
eH = eventHistory.eventHistory()
eH.outFile("Scratch/PionEventInstanceTestOut.root")
eH.rootStructure()
runNum = 149
s_pi = 0.0
z = -nuTrLnCmplx.TrfLineCmplxLen()
t = 0.
weight = 1/len(piEI)
pdg = pc.pdgCode()
for i in range(len(piEI)):
    eventNum = i
    p = nuStrt.GeneratePiMmtm(Ppion)
    x, y, xp, yp = nuStrt.GenerateTrans(s)
    ##Use the following 3 lines instead of the 2 above to create pion characteristics according to histogram distributions
    #x, xp = RndmGen.getRandom2D()
    #y, yp = RndmGen.getRandom2D()
    #p  = RndmGen.getRandom()
    pz = np.sqrt(p**2/(1+xp**2+yp**2))
    px = xp*pz
    py = yp*pz
    pz = np.sqrt(p**2-px**2-py**2)
    pi = particle.particle(runNum,eventNum,s_pi,x,y,z,px,py,pz,t,weight,pdg)
    eH.addParticle('target', pi)
    eH.fill()
    piEI_h.append(piEvtInst.PionEventInstance(particleTarLocal=pi))
for i in range(5):
    print("    piEI_h[i]:", piEI_h[i])
    print("    ppi:", piEI_h[i].getppiGen())

eH.cd()
eH.write()
eH.outFileClose()

#! Plot result of soak test with distributions from histogram:
PionEventInstanceTest = PionEventInstanceTest + 1
print()
print("PionEventInstanceTest:", PionEventInstanceTest, \
      " plots from soak test generating distributions from particle input.")

s_h        = np.array([])
pxpi_h      = np.array([])
pypi_h      = np.array([])
pzpi_h      = np.array([])
Emu_h      = np.array([])
pxmu_h      = np.array([])
pymu_h      = np.array([])
pzmu_h      = np.array([])
Enumu_h    = np.array([])
pxnumu_h      = np.array([])
pynumu_h      = np.array([])
pznumu_h      = np.array([])
phi_h      = np.array([])
costheta_h = np.array([])

for piEvt_h in piEI_h:
    s_h     = np.append(s_h,     piEvt_h.getTraceSpaceCoord()[0])
    pxpi_h = np.append(pxpi_h,   piEvt_h.getTraceSpaceCoord()[4]*piEvt_h.getppi())
    pypi_h = np.append(pypi_h,   piEvt_h.getTraceSpaceCoord()[5]*piEvt_h.getppi())
    pzpi_h = np.append(pzpi_h,   np.sqrt(piEvt_h.getppiGen()**2 - pxpi_h**2 - pypi_h**2))
    Emu_h = np.append(Emu_h, piEvt_h.getmu4mmtm()[0])
    pxmu_h = np.append(pxmu_h, piEvt_h.getmu4mmtm()[1][0])
    pymu_h = np.append(pymu_h, piEvt_h.getmu4mmtm()[1][1])
    pzmu_h = np.append(pzmu_h, piEvt_h.getmu4mmtm()[1][2])
    Enumu_h = np.append(Enumu_h, piEvt_h.getnumu4mmtm()[0])
    pxnumu_h = np.append(pxnumu_h, piEvt_h.getnumu4mmtm()[1][0])
    pynumu_h = np.append(pynumu_h, piEvt_h.getnumu4mmtm()[1][1])
    pznumu_h = np.append(pznumu_h, piEvt_h.getnumu4mmtm()[1][2])
    phi_h = np.append(phi_h, piEvt_h.getphi())
    costheta_h = np.append(costheta_h, piEvt_h.getcostheta())


dir = os.path.join(nuSIMPATH,r'Scratch/PionEventInstanceTst')
dirExist = os.path.isdir(dir)
if dirExist == False:
    os.mkdir(dir)

'''
#Plots from regular method using ppi
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
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_s.pdf')
plt.close()

#-- Energy distributions:
n, bins, patches = plt.hist(Emu, bins=50, color='y', range=(0.,10.0))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('muon energy distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_Emu.pdf')
plt.close()

n, bins, patches = plt.hist(Enumu, bins=50, color='y', range=(0.,10.))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('Muon-neutrino energy distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_Enumu.pdf')
plt.close()

#-- Momentum distributions:
n, bins, patches = plt.hist(pxpi, bins=50, color='y')
plt.xlabel('Px (GeV)')
plt.ylabel('Frequency')
plt.title('pion px distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_PxPi.pdf')
plt.close()

n, bins, patches = plt.hist(pypi, bins=50, color='y')
plt.xlabel('Py (GeV)')
plt.ylabel('Frequency')
plt.title('pion py distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_PyPi.pdf')
plt.close()

n, bins, patches = plt.hist(pzpi, bins=50, color='y', range=(4.,6.))
plt.xlabel('Pz (GeV)')
plt.ylabel('Frequency')
plt.title('pion pz distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_PzPi.pdf')
plt.close()

n, bins, patches = plt.hist(pxmu, bins=50, color='y')
plt.xlabel('Px (GeV)')
plt.ylabel('Frequency')
plt.title('muon px distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_PxMu.pdf')
plt.close()

n, bins, patches = plt.hist(pymu, bins=50, color='y')
plt.xlabel('Py (GeV)')
plt.ylabel('Frequency')
plt.title('muon py distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_PyMu.pdf')
plt.close()

n, bins, patches = plt.hist(pzmu, bins=50, color='y')
plt.xlabel('Pz (GeV)')
plt.ylabel('Frequency')
plt.title('muon pz distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_PzMu.pdf')
plt.close()

n, bins, patches = plt.hist(pxnumu, bins=50, color='y')
plt.xlabel('Px (GeV)')
plt.ylabel('Frequency')
plt.title('muon neutrino px distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_PxNumu.pdf')
plt.close()

n, bins, patches = plt.hist(pynumu, bins=50, color='y')
plt.xlabel('Py (GeV)')
plt.ylabel('Frequency')
plt.title('muon neutrino py distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_PyNumu.pdf')
plt.close()

n, bins, patches = plt.hist(pznumu, bins=50, color='y')
plt.xlabel('Pz (GeV)')
plt.ylabel('Frequency')
plt.title('muon neutrino pz distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_PzNumu.pdf')
plt.close()

n, bins, patches = plt.hist(phi, bins=50, color='y', range=(0.,6.3))
plt.xlabel('Angle')
plt.ylabel('Frequency')
plt.title('Muon phi distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_phi.pdf')
plt.close()

n, bins, patches = plt.hist(costheta, bins=50, color='y', range=(-2.,2.))
plt.xlabel('$\\cos(\\theta)$')
plt.ylabel('Frequency')
plt.title('Muon $\\cos(\\theta)$ distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst5_cosTheta.pdf')
plt.close()

#Plots from method using eventHistory
#-- Lifetime distribution:
n, bins, patches = plt.hist(s_h, bins=50, color='y', log=True)
plt.xlabel('s (m)')
plt.ylabel('Entries')
plt.title('s distribution')
## add a 'best fit' line - convert en
charDist = Ppion*1000.0*pc.SoL()*pc.lifetime()/pc.mass()
slope = 1./charDist
y = n[0]*np.exp(-slope*bins)
plt.plot(bins, y, '-', color='b')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_s.pdf')
plt.close()

#-- Energy distributions:
n, bins, patches = plt.hist(Emu_h, bins=50, color='y', range=(0.,10.0))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('muon energy distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_Emu.pdf')
plt.close()

n, bins, patches = plt.hist(Enumu_h, bins=50, color='y', range=(0.,10.))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('Muon-neutrino energy distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_Enumu.pdf')
plt.close()

#-- Momentum distributions:
n, bins, patches = plt.hist(pxpi_h, bins=50, color='y')
plt.xlabel('Px (GeV)')
plt.ylabel('Frequency')
plt.title('pion px distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_PxPi.pdf')
plt.close()

n, bins, patches = plt.hist(pypi_h, bins=50, color='y')
plt.xlabel('Py (GeV)')
plt.ylabel('Frequency')
plt.title('pion py distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_PyPi.pdf')
plt.close()

n, bins, patches = plt.hist(pzpi_h, bins=50, color='y', range=(4.,6.))
plt.xlabel('Pz (GeV)')
plt.ylabel('Frequency')
plt.title('pion pz distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_PzPi.pdf')
plt.close()

n, bins, patches = plt.hist(pxmu_h, bins=50, color='y')
plt.xlabel('Px (GeV)')
plt.ylabel('Frequency')
plt.title('muon px distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_PxMu.pdf')
plt.close()

n, bins, patches = plt.hist(pymu_h, bins=50, color='y')
plt.xlabel('Py (GeV)')
plt.ylabel('Frequency')
plt.title('muon py distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_PyMu.pdf')
plt.close()

n, bins, patches = plt.hist(pzmu_h, bins=50, color='y')
plt.xlabel('Pz (GeV)')
plt.ylabel('Frequency')
plt.title('muon pz distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_PzMu.pdf')
plt.close()

n, bins, patches = plt.hist(pxnumu_h, bins=50, color='y')
plt.xlabel('Px (GeV)')
plt.ylabel('Frequency')
plt.title('muon neutrino px distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_PxNumu.pdf')
plt.close()

n, bins, patches = plt.hist(pynumu_h, bins=50, color='y')
plt.xlabel('Py (GeV)')
plt.ylabel('Frequency')
plt.title('muon neutrino py distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_PyNumu.pdf')
plt.close()

n, bins, patches = plt.hist(pznumu_h, bins=50, color='y')
plt.xlabel('Pz (GeV)')
plt.ylabel('Frequency')
plt.title('muon neutrino pz distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_PzNumu.pdf')
plt.close()


n, bins, patches = plt.hist(phi_h, bins=50, color='y', range=(0.,6.3))
plt.xlabel('Angle')
plt.ylabel('Frequency')
plt.title('Muon phi distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_phi.pdf')
plt.close()

n, bins, patches = plt.hist(costheta_h, bins=50, color='y', range=(-2.,2.))
plt.xlabel('$\\cos(\\theta)$')
plt.ylabel('Frequency')
plt.title('Muon $\\cos(\\theta)$ distribution')
#plt.show()
plt.savefig('Scratch/PionEventInstanceTst/PionEventInstanceTst6_cosTheta.pdf')
plt.close()
'''
##! Complete:
print()
print("========  PionEventInstance: tests complete  ========")
