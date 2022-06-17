#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for going through the eventHistory and checking simulation
=================================================================

  Assumes that nuSim code is in python path.

    @file checkGeneration.py

    @brief   anaylyse the eventHistory

    @author  Paul Kyberd

    @version     1.0
    @date        06 June 2022


"""

# Generic Python imports
import sys, os, argparse

from pathlib import Path            # checking for file existance
import csv                          # so I can read a synthetic data file
import math
import numpy

# nuStorm imports
import eventHistory as eventHistory
import histoManager as histoManager
import histsCreate as histsCreate
import particle as particle
import control
import logging

# method to parse arguments
class keyvalue(argparse.Action):
  # Constructor
  def __call__( self, parser, namespace, values, option_string = None):
    d = getattr(namespace, self.dest) or {}

    if values:
            for item in values:
                split_items = item.split("=", 1)
                key = split_items[
                    0
                ].strip()  # we remove blanks around keys, as is logical
                value = split_items[1]

                d[key] = value

    setattr(namespace, self.dest, d)

##! Start:

aHVersion = 1.0;

#  if there are any keys passed in the command line put them in a dictionary
parser = argparse.ArgumentParser()
parser.add_argument("--set", metavar="KEY=VALUE", 
                    nargs='+', 
                    help="Set a number of key-VALUEW PAIRS", 
                    action=keyvalue)
args = parser.parse_args()
print(args.set)
#  set a flag only if the parameters dictions was created
runNumberFlg = False
if args.set:
  if "runNumber" in args.set.keys():
    print ("runNumber  in dict")
    runNumberFlg = True 
    runNumberVal = args.set["runNumber"]
  else:
    print ("runNumber not in dict")

nTests = 0
testFails = 0
descriptions=[]
locations=[]
testStatus=[]

# Pions Flash in the production straight
#controlFile = "102-Studies/pencilValidation/PSPiFlash.dict"
# Muons decay in ring, after first pass of the production straight
#controlFile = "102-Studies/pencilValidation/PSMuRingDcy.dict"
# Muons decay in production straight with pions.
#controlFile = "102-Studies/pencilValidation/TLPiFlash.dict"
#controlFile = "102-Studies/pencilValidation/PSMuDcy.dict"
# muons from transfer line
StudyDir = os.getenv('StudyDir')
StudyName = os.getenv('StudyName')
controlFile = os.path.join(StudyDir, StudyName,"PSPiFLash3.dict")
ctrlInst = control.control(controlFile)

#       logfile initialisation
logging.basicConfig(filename=ctrlInst.logFile(), encoding='utf-8', level=logging.INFO)
print("========  checks on the eventHistory start  ========")
logging.info("\n\n")
logging.info("========  checks on the eventHistory: start  ======== Version %s", aHVersion)
logging.info("Control file %s", controlFile)

objRd = eventHistory.eventHistory()
# Get the nuSIM path name and use it to set names for the inputfile and the outputFile
nuSIMPATH = os.getenv('nuSIMPATH')

if runNumberFlg:
  rootFilename = os.path.join(StudyDir, StudyName, 'normalisation' + str(runNumberVal)+'.root')
else:
  rootFilename = os.path.join(StudyDir, StudyName, 'normalisation' + str(ctrlInst.runNumber())+'.root')

print(f"rootFileName is {rootFilename}")

logging.info("  input file %s ", rootFilename)

objRd.inFile(rootFilename)

nEvent = objRd.getEntries()
print ("number of entries is ", nEvent)


dbgFlag = False

print("start reading")
# Now read the input from the csv file, create the eventHistory and com[are it to the one recovered from the root file
#  get the first event fromnext event from root
locRecord=[]
locStatus=[]
hm = histoManager.histoManager()
hmDataOut = histoManager.histoManager()
# common histograms at all the points in the history
hC = histsCreate.histsCreate(hm, ctrlInst.plotsDict())
hC.histAdd("target")

errors = 0
piDcyXPos=numpy.array([0,0])    # Postion flags
piDcyYPos=numpy.array([0,0])
piDcyZPos=numpy.array([0,0])    
piDcyPx=numpy.array([0,0])     # Momentum flags
piDcyPy=numpy.array([0,0])
piDcyPz=numpy.array([0,0])
piDcyMassShell=numpy.array([0,0])
piDcyMuMassShell=numpy.array([0,0])
piDcyNuMassShell=numpy.array([0,0])
piDcyE=numpy.array([0,0])

for pnt in range(nEvent):
# read an event
    if (pnt%1000 == 0): print ("Event ", pnt)
    if (pnt > 2): break
    objRd.readNext()
#    objRd.display()
    partTar = objRd.findParticle("target")
    hC.histsFill("target", partTar)
    if (dbgFlag): print ("target Particle is ", partTar)

    partPS = objRd.findParticle("productionStraight")
    if (dbgFlag): print ("productionStraight Particle is ", partPS)

    partPSEnd = objRd.findParticle("prodStraightEnd")
    if (dbgFlag): print ("prodStraightEnd Particle is ", partPSEnd)

    partPD = objRd.findParticle("pionDecay")
    if (dbgFlag): print ("pionDecay Particle is ", partPD)
    if (partPD.weight() >  0.0):
      if (dbgFlag): print (f"pion Decay {partPD}")

    partMuonP = objRd.findParticle("muonProduction")
    if (dbgFlag): print ("muonProduction Particle is ", partMuonP)
    if (partPD.weight() >  0.0):
      if (dbgFlag): print (f"muon Production {partMuonP}")

    partFlshNu = objRd.findParticle("piFlashNu")
    if (dbgFlag): print ("piFlashNu Particle is ", partFlshNu)
    if (partPD.weight() >  0.0):
      if (dbgFlag): print (f"Pion Flash neutrino {partFlshNu}")

    partMuonDecay = objRd.findParticle("muonDecay")
    if (dbgFlag): print ("muondecay Particle is ", partMuonDecay)

    partEProd = objRd.findParticle("eProduction")
    if (dbgFlag): print ("eProduction Particle is ", partEProd)

    partNumuProd = objRd.findParticle("numuProduction")
    if (dbgFlag): print ("numuProduction Particle is ", partNumuProd)

    partNueProd = objRd.findParticle("nueProduction")
    if (dbgFlag): print ("nueProduction Particle is ", partNueProd)

    partNumuDetect = objRd.findParticle("numuDetector")
    if (dbgFlag): print ("numuProduction Particle is ", partNumuDetect)

    partNueDetect = objRd.findParticle("nueDetector")
    if (dbgFlag): print ("nueProduction Particle is ", partNueDetect)

#  Check the pion decay kinematics - if there is a decay
    if partPD.weight() > 0:
#     First the decay position - can do a straight equals because they are just copied
      if (partPD.x() != partMuonP.x()) or (partPD.x() != partFlshNu.x()):
        print ("error in the x position of the pion decay")
        piDcyXPos[0] = piDcyXPos[0] + 1
        errors = errors + 1
      else:
        piDcyXPos[1] = piDcyXPos[1] + 1
      if (partPD.y() != partMuonP.y()) or (partPD.y() != partFlshNu.y()):
        print ("error in the y position of the pion decay")
        piDcyYPos[0] = piDcyYPos[0] + 1
        errors = errors + 1
      else:
        piDcyYPos[1] = piDcyYPos[1] + 1
      if (partPD.z() != partMuonP.z()) or (partPD.z() != partFlshNu.z()):
        print ("error in the y position of the pion decay")
        piDcyZPos[0] = piDcyZPos[0] + 1
        errors = errors + 1
      else:
        piDcyZPos[1] = piDcyZPos[1] + 1

#  now momentum conservation ... have to give possible error here - go with 1 in 10^6
      print(f"P pion {partPD.p()}")
      print(f"P neutrino {partFlshNu.p()}")
      print(f"P muon {partMuonP.p()}")
      if abs(partMuonP.p()[1][0] + partFlshNu.p()[1][0] - partPD.p()[1][0]) < 1e-6:
        piDcyPx[1] = piDcyPx[1] + 1
      else:
        print ("error in the Px of the pion decay")
        piDcyPx[0] = piDcyPx[0] + 1
        errors = errors + 1
      if abs(partMuonP.p()[1][1] + partFlshNu.p()[1][1] - partPD.p()[1][1]) < 1e-6:
        piDcyPy[1] = piDcyPy[1] + 1
      else:
        print ("error in the Py of the pion decay")
        piDcyPy[0] = piDcyPy[0] + 1
        errors = errors + 1
      deltaPz = partMuonP.p()[1][2] + partFlshNu.p()[1][2] - partPD.p()[1][2]
      if abs(deltaPz) < 1e-6:
        piDcyPz[1] = piDcyPz[1] + 1
      else:
        print (f"error in the Pz of the pion decay. Delta Pz is {deltaPz}")
        piDcyPz[0] = piDcyPz[0] + 1
        errors = errors + 1
#  Energy conservation at the pion decay
      deltaE = partPD.p()[0] - partMuonP.p()[0] - partFlshNu.p()[0]
      if abs(deltaE) < 1e-6:
        piDcyE[1] = piDcyE[1] + 1
      else:
        print ("error Energy not conserved at pion decay")
        piDcyE[0] = piDcyE[0] + 1
        print(f"energy mismatch {deltaE}")
        errors = errors + 1
        print(f"dE + dPz {deltaE+deltaPz}")
#  Pion on the mass shell to 1 in 10^6
      deltaMPi = partPD.p()[1][0]**2 + partPD.p()[1][1]**2 + partPD.p()[1][2]**2 + partPD.mass()**2 - partPD.p()[0]**2
      if abs(deltaMPi) < 1e-6:
        piDcyMassShell[1] = piDcyMassShell[1] + 1
      else:
        print ("error pion not on the mass shell at pion decay")
        piDcyMassShell[0] = piDcyMassShell[0] + 1
        print(f"mass deficit {deltaMPi}")
        errors = errors + 1
        exit()
#  Muon the mass shell to 1 in 10^6
      deltaMMu = partMuonP.p()[1][0]**2 + partMuonP.p()[1][1]**2 + partMuonP.p()[1][2]**2 + partMuonP.mass()**2 - partMuonP.p()[0]**2
      if abs(deltaMMu) < 1e-6:
        piDcyMuMassShell[1] = piDcyMuMassShell[1] + 1
      else:
        print ("error muon not on the mass shell at pion decay point")
        piDcyMuMassShell[0] = piDcyMuMassShell[0] + 1
        print(f"mass deficit {deltaMMu}")
        errors = errors + 1
#  neutrino the mass shell to 1 in 10^6
      deltaMNu = partFlshNu.p()[1][0]**2 + partFlshNu.p()[1][1]**2 + partFlshNu.p()[1][2]**2 + partFlshNu.mass() - partFlshNu.p()[0]**2
      if abs(deltaMNu) < 1e-6:
        piDcyNuMassShell[1] = piDcyNuMassShell[1] + 1
      else:
        print ("error neutrino not on the mass shell at pion decay point")
        piDcyNuMassShell[0] = piDcyNuMassShell[0] + 1
        print(f"mass deficit {deltaMNu}")
        errors = errors + 1

#  Check the muon decay kinematics - if there is a decay
        if partMuonDecay.weight() > 1:
          print ("found a muon decay")
          exit()
print()
print (f"piDcyXpos: failed {piDcyXPos[0]} passed  {piDcyXPos[1]} ")
print (f"piDcyYpos {piDcyYPos}")
print (f"piDcyZpos {piDcyZPos}")
print (f"piDcyPx {piDcyPx}")
print (f"piDcyPy {piDcyPy}")
print (f"piDcyPz {piDcyPz}")
print (f"Energy conservation at pion decay {piDcyE}")
print (f"pion decay pion on mass shell {piDcyMassShell}")
print (f"pion decay muon on mass shell {piDcyMuMassShell}")
print (f"pion decay neutrino on mass shell {piDcyNuMassShell}")
print (f"Total number of errors {errors}")
print("========  checking the eventHistory ends ========")
logging.info("========  checking the eventHistory: end  ========")


sys.exit(0)
