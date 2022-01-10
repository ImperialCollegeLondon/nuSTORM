#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "RandomGenerator" class ... initialisation and get methods
==================================

  RandomGenerator.py -- set "relative" path to code

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import RandomGenerator as Rndm

##! Start:
print("========  RandomGenerator: tests start  ========")

##! Test singleton class feature:
RandomGeneratorTest = 1
print()
print("RandomGeneratorTest:", RandomGeneratorTest, " check if class is a singleton.")
nuSIMPATH    = os.getenv('nuSIMPATH')
rootfilename = os.path.join(nuSIMPATH, 'Scratch/plots_endProdStrght.root')
print(rootfilename)
RndmGen  = Rndm.RandomGenerator(rootfilename)
RndmGen1 = Rndm.RandomGenerator(rootfilename)
print("---->RndmGen singleton test:", id(RndmGen), id(RndmGen1), id(RndmGen)-id(RndmGen1))
if RndmGen != RndmGen1:
    raise Exception("RandomGenerator is not a singleton class!")

##! Check get methods:
RandomGeneratorTest = 2
print()
print("RandomGeneratorTest:", RandomGeneratorTest, " check get methods.")
print("    RandomGenerator: version:", RndmGen.CdVrsn())

x = np.array([])
y = np.array([])
xp = np.array([])
yp = np.array([])
p = np.array([])

n = 200000

for i in range(0,n):
    x_cur, xp_cur = RndmGen.getRandom2D('histXPS')
    y_cur, yp_cur = RndmGen.getRandom2D('histYPS')
    x = np.append(x,x_cur)
    y = np.append(y,y_cur)
    xp = np.append(xp,xp_cur)
    yp = np.append(yp,yp_cur)
    p_cur = RndmGen.getRandom('histP')
    p = np.append(p,p_cur)

n, bins, patches = plt.hist(p, bins=50, color='y', range=(4.25,5.75))
plt.xlabel('Momentum (GeV)')
plt.ylabel('Entries')
plt.title('Random momentum distribution')
plt.savefig('Scratch/RndmGen_p_rndm.pdf')
plt.close()

pltX = plt.hist2d(x,xp, bins=200,range=[[-0.5,0.5],[-.05,.05]],norm=colors.LogNorm())
plt.colorbar()
plt.xlabel('x [m]')
plt.ylabel('x^prime')
plt.title('Random phase space distribution in x')
plt.savefig('Scratch/RndmGen_xPS_rndm.pdf')
plt.close()

pltY = plt.hist2d(y,yp, bins=200, range=[[-0.5,0.5],[-.05,.05]],norm=colors.LogNorm())
plt.colorbar()
plt.xlabel('y [m]')
plt.ylabel('y^prime')
plt.title('Random phase space distribution in y')
plt.savefig('Scratch/RndmGen_yPS_rndm.pdf')
plt.close()

##! Complete:
print()
print("========  RandomGenerator: tests complete  ========")
