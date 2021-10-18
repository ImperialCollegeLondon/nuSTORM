#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "nuSTORMTrfLine" class
========================================

  Assumes python path includes nuSim code.

  Script starts by testing built in methods.  Then a soak test with a set
  of reference plots.

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from Simulation import *
import nuSTORMTrfLine as nuTrf

##! Start:
print("========  nuSTORMTrfLine: tests start  ========")

##! Test singleton class feature:
nuSTORMTrfLineTest = 1
print()
print("nuSTORMTrfLineTest:", nuSTORMTrfLineTest, " check if class is a singleton.")
nuSIMPATH = os.getenv('nuSIMPATH')
#filename  = os.path.join(nuSIMPATH,r'11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
filename = r'/home/marvin/Documents/masters_thesis/nuSTORM/11-Parameters/nuSTORM-TrfLine-Params-v1.0.csv'
nuTrfLine  = nuTrf.nuSTORMTrfLine(filename)
nuTrfLine1 = nuTrf.nuSTORMTrfLine(filename)
print("    nuTrfLine singleton test:", id(nuTrfLine), id(nuTrfLine1), id(nuTrfLine)-id(nuTrfLine1))
if nuTrfLine != nuTrfLine1:
    raise Exception("nuSTORMTrfLine is not a singleton class!")

##! Check built-in methods:
nuSTORMTrfLineTest = 2
print()
print("nuSTORMTrfLineTest:", nuSTORMTrfLineTest, " check built-in methods.")
print("    __repr__:")
print(nuTrfLine)

##! Check get methods:
nuSTORMTrfLineTest = 3
print()
print("nuSTORMTrfLineTest:", nuSTORMTrfLineTest, " check get methods.")
print("----> printParams() method; reports parameters loaded")
nuTrfLine.printParams()

##! Check momentum, z, and transverse distributions:
nuSTORMTrfLineTest = 4
print()
print("nuSTORMTrfLineTest:", nuSTORMTrfLineTest, " momentum distribution.")
pmu = np.array([])
nmu = 0
zmu = np.array([])
xmu = np.array([])
ymu = np.array([])
xpmu = np.array([])
ypmu = np.array([])
print()
for i in range(500):   #50000
  nmu += 1
  pmu = np.append(pmu, nuTrfLine.GenerateMmtm(5.))
  zmu = np.append(zmu, nuTrfLine.Calculatez(-50.*getRandom()))
  x, y, xp, yp = nuTrfLine.GenerateTrans(-50.*getRandom())
  xmu  = np.append(xmu, x)
  ymu  = np.append(ymu, y)
  xpmu = np.append(xpmu, xp)
  ypmu = np.append(ypmu, yp)
  if nmu < 11:
      print(nmu, pmu[nmu-1], zmu[nmu-1], x, y, xp, yp)

##! Check generating pion distributions from scratch
#nuSTORMTrfLineTest = 5
#print()
#print("nuSTORMTrfLineTest:", nuSTORMTrfLineTest, " pion distributions.")
#pis = []
#Epi = np.array([])
#ppi = np.array([])
#xppi = np.array([])
#yppi = np.array([])
#sd = np.array([])
#p0 = 5.
#s_final=-50.
#runNum = 101
#eventNum = 0
#print()
#for i in range(500000):
#    eventNum +=1
#    s_decay, pi = nuTrfLine.GeneratePion(p0=p0,s_final=s_final,runNum=runNum,eventNum=eventNum)
#    pis.append(pi)
#    Epi = np.append(Epi,pi.p()[0])
#    sd = np.append(sd,s_decay)
#    if eventNum < 11:
#        print(eventNum,pi.run(),s_decay,pi.s(),pi.x(),pi.y(),pi.z(),pi.p()[1][0],pi.p()[1][1],pi.p()[1][2],pi.mass(),pi.t(),pi.weight())
#    elif eventNum%100000 == 0:
#        print(eventNum)


##! PLOTTING
#For nuSTORMTrfLineTest = 4
n, bins, patches = plt.hist(pmu, bins=50, color='y', range=(4.0,6.0))
plt.xlabel('Momentum (GeV)')
plt.ylabel('Entries')
plt.title('momentum distribution')
plt.savefig('Scratch/nuSTORMTrfLineTst_plot1.pdf')
plt.close()

n, bins, patches = plt.hist(zmu, bins=50, color='y', range=(-50.,0.))
plt.xlabel('z (m)')
plt.ylabel('Entries')
plt.title('z distribution')
plt.savefig('Scratch/nuSTORMTrfLineTst_plot2.pdf')
plt.close()

n, bins, patches = plt.hist(xmu, bins=50, color='y', range=(-0.15,0.15))
plt.xlabel('x (m)')
plt.ylabel('Entries')
plt.title('x distribution')
plt.savefig('Scratch/nuSTORMTrfLineTst_plot3.pdf')
plt.close()

n, bins, patches = plt.hist(ymu, bins=50, color='y', range=(-0.15,0.15))
plt.xlabel('y (m)')
plt.ylabel('Entries')
plt.title('y distribution')
plt.savefig('Scratch/nuSTORMTrfLineTst_plot4.pdf')
plt.close()

n, bins, patches = plt.hist(xpmu, bins=50, color='y', range=(-0.0075,0.0075))
plt.xlabel('x^prime')
plt.ylabel('Entries')
plt.title('x^prime distribution')
plt.savefig('Scratch/nuSTORMTrfLineTst_plot5.pdf')
plt.close()

n, bins, patches = plt.hist(ypmu, bins=50, color='y', range=(-0.0075,0.0075))
plt.xlabel('y^prime')
plt.ylabel('Entries')
plt.title('y^prime distribution')
plt.savefig('Scratch/nuSTORMTrfLineTst_plot6.pdf')
plt.close()

#For nuSTORMTrfLineTest = 5
n, bins, patches = plt.hist(Epi, bins=50, color='y', range=(4.0,6.0))
plt.xlabel('Energy (GeV)')
plt.ylabel('Entries')
plt.title('energy distribution')
plt.savefig('Scratch/nuSTORMTrfLineTst_pi_plot1.pdf')
plt.close()

n, bins, patches = plt.hist(sd, bins=50, color='y', range=(0.,10e-9))
plt.xlabel('z decay (m)')
plt.ylabel('Entries')
plt.title('z decay distribution')
plt.savefig('Scratch/nuSTORMTrfLineTst_pi_plot2.pdf')
plt.close()

##! Complete:
print()
print("========  nuSTORMTrfLine: tests complete  ========")
