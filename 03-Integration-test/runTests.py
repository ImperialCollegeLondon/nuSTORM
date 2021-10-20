#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class testDrive:
=============

runs the test scripts and picks up the return codes

V1.2 							20 October 2021 						
	Add RunSimulation.py

V1.1 							20 October 2021 						
	Add test for PionEventInstance.py

Version 1.0 				A framework and some of the tests



"""
from subprocess import run

tests=[]
notes=[]
testDesc=[]
testStatus=[]
testPnt = 0
testFails = 0

#		Test 1
# test muonConst
testFile = '02-Tests/MuonConstTst.py'
testDesc.append(testFile)
notes.append("Need to check output by hand to complete tests")
test = run(['python', testFile])
# This returns a 1 if it the class is not a singleton - else need to check the 
# output
if (test.returncode == 0):
	print (testFile, " passes")
	print (notes[testPnt])
	testStatus.append("passes")
else:
	print (testFile, " fails")
	testFails = testFails + 1
	testStatus.append("fails")

testPnt = testPnt + 1
#		Test 2
# test pionConst
testFile = '02-Tests/PionConstTst.py'
testDesc.append(testFile)
notes.append("for fail check Scratch/pionTst.out against 02-Tests/referenceOutput/pionTst.ref")
test = run(['python', testFile])
# This returns a 1 if it the class is not a singleton - else need to check the 
# output
if (test.returncode == 0):
	print (testFile, " passes")
	print (notes[testPnt])
	testStatus.append("passes")
else:
	print (testFile, " fails")
	testFails = testFails + 1
	testStatus.append("fails")

testPnt = testPnt + 1
#		Test 3
# test traceSpace
testFile = '02-Tests/traceSpaceTst.py'
testDesc.append(testFile)
notes.append("__stre__ and __repr__ - need checking in the print out")
test = run(['python', testFile])

if (test.returncode == 0):
	print (testFile, " passes")
	print (notes[testPnt])
	testStatus.append("passes")

else:
	print (testFile, " fails")
	testFails = testFails + 1
	testStatus.append("fails")

testPnt = testPnt + 1
#		Test 4
# test particle
testFile = '02-Tests/particleTst.py'
testDesc.append(testFile)
notes.append("for fail check Scratch/particleTst.out against 02-Tests/particleTst.ref")
test = run(['python', testFile])

if (test.returncode == 0):
	print (testFile, " passes")
	print (notes[testPnt])
	testStatus.append("passes")

else:
	print (testFile, " fails")
	testFails = testFails + 1
	testStatus.append("fails")

testPnt = testPnt + 1
#		Test 5
# test pion
testFile = '02-Tests/pionTst.py'
testDesc.append(testFile)
notes.append("On fail check Scratch/particleTst.out against 02-Tests/particleTst.ref")
test = run(['python', testFile])

if (test.returncode == 0):
	print (testFile, " passes")
	print (notes[testPnt])
	testStatus.append("passes")

else:
	print (testFile, " fails")
	testFails = testFails + 1
	testStatus.append("fails")

testPnt = testPnt + 1
#		Test 5
# test muon decay
testFile = '02-Tests/MuonDecayTst.py'
testDesc.append(testFile)
notes.append("not run yet")
#test = run(['python', testFile])

#if (test.returncode == 0):
#	print (testFile, " passes")
#	print (notes[testPnt])
#	testStatus.append("passes")
#
#else:
#	print (testFile, " fails")
#	testFails = testFails + 1
#	testStatus.append("fails")
testStatus.append("Not run")

testPnt = testPnt + 1
#		Test 6
# test truncatedmuon decay
testFile = '02-Tests/MuonDecayTrncTst.py'
testDesc.append(testFile)
notes.append("not run yet")
#test = run(['python', testFile])

#if (test.returncode == 0):
#	print (testFile, " passes")
#	print (notes[testPnt])
#	testStatus.append("passes")
#
#else:
#	print (testFile, " fails")
#	testFails = testFails + 1
#	testStatus.append("fails")
testStatus.append("Not run")

testPnt = testPnt + 1
#		Test 7
# test truncatedmuon decay
testFile = '02-Tests/nuSTORMPrdStrghtTst.py'
testDesc.append(testFile)
notes.append("not run yet")
#test = run(['python', testFile])

#if (test.returncode == 0):
#	print (testFile, " passes")
#	print (notes[testPnt])
#	testStatus.append("passes")
#
#else:
#	print (testFile, " fails")
#	testFails = testFails + 1
#	testStatus.append("fails")
testStatus.append("Not run")

testPnt = testPnt + 1
#		Test 8
# test truncatedmuon decay
testFile = '02-Tests/NeutrinoEventInstanceTst.py'
testDesc.append(testFile)
notes.append("not run yet")
#test = run(['python', testFile])

#if (test.returncode == 0):
#	print (testFile, " passes")
#	print (notes[testPnt])
#	testStatus.append("passes")
#
#else:
#	print (testFile, " fails")
#	testFails = testFails + 1
#	testStatus.append("fails")
testStatus.append("Not run")

testPnt = testPnt + 1
#		Test 9
# test truncatedmuon decay
testFile = '02-Tests/MuonBeam4CoolingDemoTst.py'
testDesc.append(testFile)
notes.append("not run yet")
#test = run(['python', testFile])

#if (test.returncode == 0):
#	print (testFile, " passes")
#	print (notes[testPnt])
#	testStatus.append("passes")
#
#else:
#	print (testFile, " fails")
#	testFails = testFails + 1
#	testStatus.append("fails")
testStatus.append("Not run")

testPnt = testPnt + 1
#		Test 10
# test truncatedmuon decay
testFile = '02-Tests/SimulationTst.py'
testDesc.append(testFile)
notes.append("Checks values of seed and random, not distribution of values")
test = run(['python', testFile])

if (test.returncode == 0):
	print (testFile, " passes")
	print (notes[testPnt])
	testStatus.append("passes")
else:
	print (testFile, " fails")
	testFails = testFails + 1
	testStatus.append("fails")

testPnt = testPnt + 1
#		Test 11
# test truncatedmuon decay
testFile = '02-Tests/RunSimulation.py'
testDesc.append(testFile)
notes.append("Just runs RunSimulation ... no checks on output - need a RunSimulationTst.py")
test = run(['python', testFile])

if (test.returncode == 0):
	print (testFile, " passes")
	print (notes[testPnt])
	testStatus.append("passes")

else:
	print (testFile, " fails")
	testFails = testFails + 1
	testStatus.append("fails")


testPnt = testPnt + 1
#		Test 12
# test truncatedmuon decay
testFile = '02-Tests/eventHistoryTst.py'
testDesc.append(testFile)
notes.append("Not checked that names eg pi+, give same results as pdg Code")
test = run(['python', testFile])

if (test.returncode == 0):
	print (testFile, " passes")
	print (notes[testPnt])
	testStatus.append("passes")

else:
	print (testFile, " fails")
	testFails = testFails + 1
	testStatus.append("fails")

testPnt = testPnt + 1
#		Test 13
# test truncatedmuon decay
testFile = '02-Tests/PionDecayTst.py'
testDesc.append(testFile)
notes.append("not run yet")
#test = run(['python', testFile])

#if (test.returncode == 0):
#	print (testFile, " passes")
#	print (notes[testPnt])
#	testStatus.append("passes")
#
#else:
#	print (testFile, " fails")
#	testFails = testFails + 1
#	testStatus.append("fails")
testStatus.append("Not run")

testPnt = testPnt + 1
#		Test 14
# test truncatedmuon decay
testFile = '02-Tests/PionDecaytrncTst.py'
testDesc.append(testFile)
notes.append("not run yet")
#test = run(['python', testFile])

#if (test.returncode == 0):
#	print (testFile, " passes")
#	print (notes[testPnt])
#	testStatus.append("passes")
#
#else:
#	print (testFile, " fails")
#	testFails = testFails + 1
#	testStatus.append("fails")
testStatus.append("Not run")

testPnt = testPnt + 1
#		Test 15
# test muon
testFile = '02-Tests/muonTst.py'
testDesc.append(testFile)
notes.append("")
test = run(['python', testFile])

if (test.returncode == 0):
	print (testFile, " passes")
	print (notes[testPnt])
	testStatus.append("passes")

else:
	print (testFile, " fails")
	testFails = testFails + 1
	testStatus.append("fails")

testPnt = testPnt + 1
#		Test 16
# test truncatedmuon decay
testFile = '02-Tests/PionDecayTst.py'
testDesc.append(testFile)
notes.append("not run yet")
#test = run(['python', testFile])

#if (test.returncode == 0):
#	print (testFile, " passes")
#	print (notes[testPnt])
#	testStatus.append("passes")
#
#else:
#	print (testFile, " fails")
#	testFails = testFails + 1
#	testStatus.append("fails")
testStatus.append("Not run")


testPnt = testPnt + 1
#		Test 17
# test truncatedmuon decay
testFile = '02-Tests/FluxCalcOutline.py'
testDesc.append(testFile)
notes.append("not run yet")
#test = run(['python', testFile])

#if (test.returncode == 0):
#	print (testFile, " passes")
#	print (notes[testPnt])
#	testStatus.append("passes")
#
#else:
#	print (testFile, " fails")
#	testFails = testFails + 1
#	testStatus.append("fails")
testStatus.append("Not run")

testPnt = testPnt + 1
#		Test 18
# test truncatedmuon decay
testFile = '02-Tests/PionEventInstanceTst.py'
testDesc.append(testFile)
notes.append("Not all output is checked automatically")
test = run(['python', testFile])

if (test.returncode == 0):
	print (testFile, " passes")
	print (notes[testPnt])
	testStatus.append("passes")

else:
	print (testFile, " fails")
	testFails = testFails + 1
	testStatus.append("fails")




print ("\n\n")
print("             SUMMARY \n")
for pnt in range(testPnt+1):
	print ("Test ", pnt, "  ",testDesc[pnt], "   ",testStatus[pnt], "             ", notes[pnt])

print ("\nTests run ", testPnt+1, " of which ", testFails," failed\n\n")

