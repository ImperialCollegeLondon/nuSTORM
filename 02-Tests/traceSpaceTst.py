#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "traceSpace" class
=================================

  Assumes that nuSim code is in python path.

  Script tests the instantiation and get methods of the class

"""

import sys
import numpy as np
import math as mt
import traceSpace as trSp

##! Start:
print("========  traceSpace: tests start  ========")

##! Create instance, test built-in methods:
traceSpaceTest = 1
print()
print("traceSpaceTest:", traceSpaceTest, " Create traceSpace and print quantities.")
sVal = 4.0
xVal = 1.1
yVal = 2.2
zVal = 9.2
xpVal = 0.44
ypVal = 0.55
tS = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal, ypVal)
print("    __str__:", tS)
print("    --repr__", repr(tS))
del tS

##! Create instance, test dynamic methods:
traceSpaceTest = 2
print()
print("traceSpaceTest:", traceSpaceTest, " Check get methods.")
tS = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal, ypVal)
if (tS.s() != sVal):
    sys.exit("traceSpace.s() is not working")
elif (tS.x() != xVal):
    sys.exit("traceSpace.x() is not working")
elif (tS.y() != yVal):
    sys.exit("traceSpace.y() is not working")
elif (tS.z() != zVal):
    sys.exit("traceSpace.z() is not working")
elif (tS.xp() != xpVal):
    sys.exit("traceSpace.xp() is not working")
elif (tS.yp() != ypVal):
    sys.exit("traceSpace.yp() is not working")
else:
    del tS

##! Complete:
print()
print("========  traceSpace:tests complete sucessfully  ========")
sys.exit(0)
