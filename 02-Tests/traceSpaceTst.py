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

testError = 0
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

tS = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal, ypVal)
tS1 = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal, ypVal)

##! Create instance, test dynamic methods:
traceSpaceTest = traceSpaceTest + 1
print("\ntraceSpaceTest:", traceSpaceTest, " check __eq__ and __ne__")
tS = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal, ypVal)
tS1 = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal, ypVal)
# check equality
if tS == tS1:
    pass
else:
    print ("traceSpaceTest: __eq__ fails")
    testError = testError + 1
del tS1

tS1 = trSp.traceSpace(sVal+0.1, xVal, yVal, zVal, xpVal, ypVal)
if tS == tS1:
    print ("traceSpaceTest: __eq__ fails")
    testError = testError + 1
del tS1

tS1 = trSp.traceSpace(sVal, xVal+0.1, yVal, zVal, xpVal, ypVal)
if tS == tS1:
    print ("traceSpaceTest: __eq__ fails")
    testError = testError + 1
del tS1

tS1 = trSp.traceSpace(sVal, xVal, yVal+0.1, zVal, xpVal, ypVal)
if tS == tS1:
    print ("traceSpaceTest: __eq__ fails")
    testError = testError + 1
del tS1

tS1 = trSp.traceSpace(sVal, xVal, yVal, zVal+0.1, xpVal, ypVal)
if tS == tS1:
    print ("traceSpaceTest: __eq__ fails")
    testError = testError + 1
del tS1

tS1 = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal+0.1, ypVal)
if tS == tS1:
    print ("traceSpaceTest: __eq__ fails")
    testError = testError + 1
del tS1

tS1 = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal, ypVal+0.1)
if tS == tS1:
    print ("traceSpaceTest: __eq__ fails")
    testError = testError + 1
del tS1

tS1 = trSp.traceSpace(sVal+0.1, xVal, yVal, zVal, xpVal, ypVal)
if tS != tS1:
    pass
else:
    print ("traceSpaceTest: __ne__ fails")
    testError = testError + 1
del tS1

tS1 = trSp.traceSpace(sVal, xVal+0.1, yVal, zVal, xpVal, ypVal)
if tS != tS1:
    pass
else:
    print ("traceSpaceTest: __ne__ fails")
    testError = testError + 1
del tS1


tS1 = trSp.traceSpace(sVal, xVal, yVal+0.1, zVal, xpVal, ypVal)
if tS != tS1:
    pass
else:
    print ("traceSpaceTest: __ne__ fails")
    testError = testError + 1
del tS1

tS1 = trSp.traceSpace(sVal, xVal, yVal, zVal+0.1, xpVal, ypVal)
if tS != tS1:
    pass
else:
    print ("traceSpaceTest: __ne__ fails")
    testError = testError + 1
del tS1

tS1 = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal+0.1, ypVal)
if tS != tS1:
    pass
else:
    print ("traceSpaceTest: __ne__ fails")
    testError = testError + 1
del tS1

tS1 = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal, ypVal+0.1)
if tS != tS1:
    pass
else:
    print ("traceSpaceTest: __ne__ fails")
    testError = testError + 1
del tS1

traceSpaceTest = traceSpaceTest + 1

print()
print("traceSpaceTest:", traceSpaceTest, " Check get methods.")
tS = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal, ypVal)
if (tS.s() != sVal):
    print("traceSpace.s() is not working")
    testError = testError + 1
elif (tS.x() != xVal):
    print("traceSpace.x() is not working")
    testError = testError + 1
elif (tS.y() != yVal):
    print("traceSpace.y() is not working")
    testError = testError + 1
elif (tS.z() != zVal):
    print("traceSpace.z() is not working")
    testError = testError + 1
elif (tS.xp() != xpVal):
    print("traceSpace.xp() is not working")
    testError = testError + 1
elif (tS.yp() != ypVal):
    print("traceSpace.yp() is not working")
    testError = testError + 1
else:
    del tS

##! Create instance and check equals
traceSpaceTest = 3
print()
print("traceSpaceTest:", traceSpaceTest, " Check equals and not equals.")
tS = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal, ypVal)
tS1 = trSp.traceSpace(sVal, xVal, yVal, zVal, xpVal, ypVal)

if (tS == tS1):
    pass
else:
    print("equality fails.")

tS2 = trSp.traceSpace(sVal+1.0, xVal, yVal, zVal, xpVal, ypVal)

if (tS != tS2):
    pass
else:
    print("not equals fails.")

##! Complete:
print()
print("========  traceSpace:tests complete   ========")
print ("\n >>>>> number of test errors is ", testError, " <<<<<<\n")
sys.exit(0)
