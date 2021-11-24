#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "PionTimeDistribution" class
=================================

  Assumes that nuSim code is in python path.

  Script starts by testing built in methods.  Then a soak test creating a
  large number of timestamps is executed.  Finally a set of reference plots
  are generated.

"""

import numpy as np
import matplotlib.pyplot as plt
import math as math
import Simulation as Simu
import PionTimeDistribution as PTD


##! Start:
print("========  PionTimeDistribution: tests start  ========")

t_bunch = 2.
t_spacing = 5.
t_extraction = 0.15

##! Create instance, test built-in methods:
PionTimeDistributionTst = 1
print()
print("PionTimeDistributionTest:", PionTimeDistributionTst, " Create pion time distribution, print quantities.")
ptd=PTD.PionTimeDistribution(t_bunch,t_spacing,t_extraction)
print("    __str__:", ptd)
print("    --repr__", repr(ptd))
del ptd

##! Create instance, test dynamic methods:
PionTimeDistributionTst = 2
print()
print("PionTimeDistributionTest:", PionTimeDistributionTst, " Create pion time distribution, print one pion timestamp.")
ptd=PTD.PionTimeDistribution(t_bunch,t_spacing,t_extraction)
print("     _t:", ptd.GenerateTime())
del ptd

##! Soak test, generate many decays:
PionTimeDistributionTst = 3
print()
print("PionTimeDistributionTest:", PionTimeDistributionTst, " Create many pion timestamps.")
ptd = PTD.PionTimeDistribution(t_bunch,t_spacing,t_extraction)
t = []
for i in range(1000000):
    t.append(ptd.GenerateTime())
for i in range(5):
    print("t = ",t[i])
    mod = math.fmod(t[i],(t_bunch+t_spacing)*10**(-9))
    print("t_mod = ",mod)
    print()

print("Plotting distributions...")
#-- Time distribution:
n, bins, patches = plt.hist(t, bins=100, color='y',range=(0.,10.*10**(-9)))
plt.xlabel('Time (s)')
plt.ylabel('Entries')
plt.title('Pion time distribution I')
plt.savefig('Scratch/PionTimeDistribution1.pdf')
plt.close()

n, bins, patches = plt.hist(t, bins=90, color='y',range=(0.,30.*10**(-9)))
plt.xlabel('Time (s)')
plt.ylabel('Entries')
plt.title('Pion time distribution II')
plt.savefig('Scratch/PionTimeDistribution2.pdf')
plt.close()

n, bins, patches = plt.hist(t, bins=150, color='y',range=(0.,t_extraction*10**(-6)))
plt.xlabel('Time (s)')
plt.ylabel('Entries')
plt.title('Pion time distribution III')
plt.savefig('Scratch/PionTimeDistribution3.pdf')
plt.close()

##! Complete:
print()
print("========  PionTimeDistribution: tests complete  ========")
