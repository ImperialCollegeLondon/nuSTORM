#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for readFluka class
=============================

  Assumes that nuSim code is in python path.


Version history:
----------------------------------------------
 1.0: 017Sept24: Test the class


"""

import sys
import os
import numpy as np
import readFluka
import histoManager

##! Start:

nTests = 0
testFails = 0
descriptions=[]
testStatus=[]
testTitle = "readFluka"

print(f"========  {testTitle} : tests start  ========")


##! Create instance, test built-in methods:
##! Create particle and print out #############################################################################
nTests = nTests + 1
descString = "Create readFluka instance"
descriptions.append(descString)

meanP = 2.0
rf = readFluka.readFluka(meanP)
#   print ("main: rf is ", rf)
print ("main: rf.curEvent is ", rf._curEvent)

rf.readEvents()
print ("number of events found ", rf.numberOfEvents())

#   print ("main: rf.inFile is ", rf._inFile)
#   print ("main: rf.rootFile is ", rf._rootFile)
#   print ("main: rf.ntuple is ", rf._ntuple)
#   print ("main: number of entries is ", rf._entries)

hm = histoManager.histoManager()
hTitle = "x"
hBins = 100
hLower = -10.0
hUpper = 10.0
xDist = hm.book(hTitle, hBins, hLower, hUpper)
print("xDist is ", xDist)

hTitle = "y"
hBins = 100
hLower = -10.0
hUpper = 10.0
yDist = hm.book(hTitle, hBins, hLower, hUpper)

hTitle = "px"
hBins = 100
hLower = -1.0
hUpper = 1.0
pxDist = hm.book(hTitle, hBins, hLower, hUpper)

hTitle = "py"
hBins = 100
hLower = -1.0
hUpper = 1.0
pyDist = hm.book(hTitle, hBins, hLower, hUpper)

hTitle = "pz"
hBins = 100
hLower = 1.7
hUpper = 2.3
pzDist = hm.book(hTitle, hBins, hLower, hUpper)

hTitle = "p"
hBins = 100
hLower = 1.7
hUpper = 2.3
pDist = hm.book(hTitle, hBins, hLower, hUpper)

hTitle = "x v y"
hLower1 = -10.0
hUpper1 = 10.0
hLower2 = -10.0
hUpper2 = 10.0
xyDist = hm.book2(hTitle, hBins, hLower1, hUpper1, hBins, hLower2, hUpper2)

hTitle = "px v py"
hLower1 = -0.6
hUpper1 = 0.6
hLower2 = -0.6
hUpper2 = 0.6
pxpyDist = hm.book2(hTitle, hBins, hLower1, hUpper1, hBins, hLower2, hUpper2)

hTitle = "x v xp"
hLower1 = -10.0
hUpper1 = 10.0
hLower2 = -0.2
hUpper2 = 0.2
xxpDist = hm.book2(hTitle, hBins, hLower1, hUpper1, hBins, hLower2, hUpper2)

hTitle = "y v yp"
hLower1 = -10.0
hUpper1 = 10.0
hLower2 = -0.2
hUpper2 = 0.2
yypDist = hm.book2(hTitle, hBins, hLower1, hUpper1, hBins, hLower2, hUpper2)

nEntries = rf.numberOfEvents()
for i in range(nEntries):
    values = rf.getNextEvent()
#   Calculate p
    px = values[2]
    py = values[3]
    pz = values[4]
    p = np.sqrt(px**2+py**2+pz**2)

    x = values[0]
    y = values[1]
    xp = px/pz
    yp = py/pz
    xDist.Fill(x)
    yDist.Fill(y)
    xyDist.Fill(x,y)
    pxDist.Fill(px)
    pyDist.Fill(py)
    pxpyDist.Fill(px,py)
    xxpDist.Fill(x, xp)
    yypDist.Fill(y, yp)
    pzDist.Fill(pz)
    pDist.Fill(p)

    if (i < 10):
        print ("event number is ", i, "  p is ", p)
        print ("rf.curEvent is ", rf._curEvent)
    elif ((i <100) and (i%10 ==0)):
        print ("event number is ", i, "  p is ", p)
        print ("rf.curEvent is ", rf._curEvent)
    elif ((i <1000) and (i%100 ==0)):
        print ("event number is ", i, "  p is ", p)
        print ("rf.curEvent is ", rf._curEvent)
    elif ((i <10000) and (i%1000 ==0)):
        print ("event number is ", i, "  p is ", p)
        print ("rf.curEvent is ", rf._curEvent)
    elif ((i <100000) and (i%10000 ==0)):
        print ("event number is ", i, "  p is ", p)
        print ("rf.curEvent is ", rf._curEvent)
    elif ((i <1000000) and (i%100000 ==0)):
        print ("event number is ", i, "  p is ", p)
        print ("rf.curEvent is ", rf._curEvent)
    elif (i%1000000 ==0):
        print ("event number is ", i, "  p is ", p)
        print ("rf.curEvent is ", rf._curEvent)



fileName="flukaReadTests.root"
hm.histOutRoot(fileName)
print(" tests finished - check file flukaReadTests.root")

##! Complete:
print()
print(f"========  {testTitle}:tests complete  ========")
print (f"\nNumber of tests is {nTests} number of fails is {testFails}")
if testFails == 0:
    sys.exit(0)
else:
    sys.exit(1)