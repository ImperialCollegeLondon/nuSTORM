

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for going through the eventHistory and making plots
==========================================================

  Assumes that nuSim code is in python path.

    @file analyseHistory.py

    @brief   anaylyse the eventHistory

    @author  Paul Kyberd

    @version     1.0
    @date        12 November 2021

    @version     1.1
    @date        19 May 2023

    Add a parser to the interface to explicitly set the run number

    @version     1.2
    @date        18 October 2024

    Remove override of nEvent ... if we reinstate make sure that the override is
    less than or equal to the number of events in the data set

    @version     1.3
    @date        18 October 2024

    Add comments to help the make sure the filling is done in the correct order



"""

# Generic Python imports
import sys, os

from pathlib import Path            # checking for file existance
import csv                          # so I can read a synthetic data file
import math
from datetime import datetime

# nuStorm imports
import eventHistory as eventHistory
import histoManager as histoManager
import histsCreate as histsCreate
import particle as particle
import control
import logging
import argparse

##! Start:

aHVersion = 1.3;

#   Define parser to allow the user to choose the run number to analyse

parser = argparse.ArgumentParser()
parser.add_argument('--run', help='Run number to use. If no run number specified the current run number is used', default='0')
args = parser.parse_args()


nTests = 0
testFails = 0
descriptions=[]
locations=[]
testStatus=[]

# muons from transfer line
StudyDir = os.getenv('StudyDir')
StudyName = os.getenv('StudyName')

# get the current run number
if (int(args.run) == 0):
  rNFile = StudyDir  + "/runNumber"
  rN = open(rNFile, "r")
  runNumber = int(rN.readline())
  rN.close()
else:
  runNumber = int(args.run)

print ("run number is ", runNumber)

# find and open logfile for reading
logFileName =os.environ['StudyDir'] + "/" + os.environ["StudyName"] + "/" + "normalisation"  + str(runNumber) + ".log"
lFile = open(logFileName, "r")
# read to get the control file used on the run that produced this log file
for line in lFile:
  index = line.find("Control File")
  if index > 0:
    start = line.find('/')
    end = line.find('.dict')
    print ("end is ---",end)
    controlFile = line[start: end+5]
    break

## control file shared with the normalisation code
ctrlInst = control.control(controlFile)


#       logfile initialisation
logging.basicConfig(filename=ctrlInst.logFile(), encoding='utf-8', level=logging.INFO)
print("========  analysing the eventHistory start  ========")
logging.info("\n\n")
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
logging.info("========  analysing the eventHistory: start  ======== Version %s... %s", aHVersion, dt_string)
logging.info("Control file %s", controlFile)

# Get the nuSIM path name and use it to set names for the inputfile and the outputFile
nuSIMPATH = os.getenv('nuSIMPATH')
rootFilename = os.path.join(StudyDir, StudyName, 'normalisation' + str(ctrlInst.runNumber())+'.root')
logging.info("  input file %s ", rootFilename)

objRd = eventHistory.eventHistory()
objRd.inFile(rootFilename)

nEvent = objRd.getEntries()
print ("number of entries is ", nEvent)
logging.info("  number of entries is  %i ", nEvent)

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
hC.histAdd("productionStraight")
hC.histAdd("prodStraightEnd")
hC.histAdd("pionDecay")
hC.histAdd("muonProduction")
hC.histAdd("piFlashNu")
hC.histAdd("muonDecay")
hC.histAdd("eProduction")
hC.histAdd("numuProduction")
hC.histAdd("nueProduction")
hC.histAdd("numuDetector")
hC.histAdd("nueDetector")
# some histograms which are not common to all points
hTitle = "flightDistNumu"
hBins = 100
hLower = 45.0
hUpper = 300.0
numuDist = hm.book(hTitle, hBins, hLower, hUpper)
hTitle = "flightDistNue"
nueDist = hm.book(hTitle, hBins, hLower, hUpper)
hTitle = "flightDistFlashNumu"
flashNumuDist = hm.book(hTitle, hBins, hLower, hUpper)
hTitle = "flightNuSourceX"
hLower = -1.0
hUpper = 1.0
flashNumuSrcX = hm.book(hTitle, hBins, hLower, hUpper)
hTitle = "flightNuSourceY"
hLower = -1.0
hUpper = 1.0
flashNumuSrcY = hm.book(hTitle, hBins, hLower, hUpper)
hTitle = "flightNuSourceZ"
hLower = -60.0
hUpper = 10.0
flashNumuSrcZ = hm.book(hTitle, hBins, hLower, hUpper)

hTitle = "time"
hLower = 0.0
hUpper = 300.0
arrivalT = hm.book(hTitle, hBins, hLower, hUpper)

hTitle = "fltNuSrcPX"
hLower = -1.0
hUpper = 1.0
flshNuSrcPX = hm.book(hTitle, hBins, hLower, hUpper)
hTitle = "fltNuSrcPY"
hLower = -1.0
hUpper = 1.0
flshNuSrcPY = hm.book(hTitle, hBins, hLower, hUpper)
hTitle = "fltNuSrcPZ"
hLower = 1.0
hUpper = 6.0
flshNuSrcPZ = hm.book(hTitle, hBins, hLower, hUpper)

strXXP = hm.book2("Prod strght: x v xp", 100, -0.1, 0.1, 100, -0.05, 0.05)
strYYP = hm.book2("Prod strght: y v yp", 100, -0.1, 0.1, 100, -0.05, 0.05)

#  Data file for genie analysis
hTitle = "Enumu"
hLower = 0.0
hUpper = 6.0
eNumuData = hmDataOut.book(hTitle, hBins, hLower, hUpper)

hTitle = "Enue"
hLower = 0.0
hUpper = 6.0
eNueData = hmDataOut.book(hTitle, hBins, hLower, hUpper)

for pnt in range(nEvent):
# read an event
    if (pnt%1000 == 0): print ("Event ", pnt)
#    if (pnt > 100): break
    objRd.readNext()
#    objRd.display()
    partTar = objRd.findParticle("target")
    hC.histsFill("target", partTar)
    if (dbgFlag): print ("target Particle is ", partTar)

    partPS = objRd.findParticle("productionStraight")
    hC.histsFill("productionStraight", partPS)
    if (dbgFlag): print ("productionStraight Particle is ", partPS)

    partPSEnd = objRd.findParticle("prodStraightEnd")
    hC.histsFill("prodStraightEnd", partPSEnd)
    if (dbgFlag): print ("prodStraightEnd Particle is ", partPSEnd)

    partPD = objRd.findParticle("pionDecay")
    hC.histsFill("pionDecay", partPD)
    if (dbgFlag): print ("pionDecay Particle is ", partPD)

    partMuonP = objRd.findParticle("muonProduction")
    hC.histsFill("muonProduction", partMuonP)
    if (dbgFlag): print ("muonProduction Particle is ", partMuonP)

    partFlshNu = objRd.findParticle("piFlashNu")
    hC.histsFill("piFlashNu", partFlshNu)
    if (dbgFlag): print ("piFlashNu Particle is ", partFlshNu)

    partMuonDecay = objRd.findParticle("muonDecay")
    hC.histsFill("muonDecay", partMuonDecay)
    if (dbgFlag): print ("muondecay Particle is ", partMuonDecay)

    partEProd = objRd.findParticle("eProduction")
    hC.histsFill("eProduction", partEProd)
    if (dbgFlag): print ("eProduction Particle is ", partEProd)

    partNumuProd = objRd.findParticle("numuProduction")
    hC.histsFill("numuProduction", partNumuProd)
    if (dbgFlag): print ("numuProduction Particle is ", partNumuProd)

    partNueProd = objRd.findParticle("nueProduction")
    hC.histsFill("nueProduction", partNueProd)
    if (dbgFlag): print ("nueProduction Particle is ", partNueProd)

    partNumuDetect = objRd.findParticle("numuDetector")
    hC.histsFill("numuDetector", partNumuDetect)
    if (dbgFlag): print ("numuProduction Particle is ", partNumuDetect)

    partNueDetect = objRd.findParticle("nueDetector")
    hC.histsFill("nueDetector", partNueDetect)
    if (dbgFlag): print ("nueProduction Particle is ", partNueDetect)

# trying hack round a problem
    sx = partPS.x()
    sxp = partPS.p()[1][0]/partPS.p()[1][2]
    strXXP.Fill(sx, sxp)
    sy = partPS.y()
    syp = partPS.p()[1][1]/partPS.p()[1][2]
    strYYP.Fill(sy, syp)

# some specfic analysis - neutrinos from muon decay if only muon decays in the ring are included
    if (partNumuDetect.weight() > 0.0):
      xs = partNumuProd.x()
      ys = partNumuProd.y()
      zs = partNumuProd.z()
      xd = partNumuDetect.x()
      yd = partNumuDetect.y()
      zd = partNumuDetect.z()
      dist = math.sqrt((xd-xs)*(xd-xs) + (yd-ys)*(yd-ys) +(zd-zs)*(zd-zs))
      numuDist.Fill(dist)
      if (abs(xd) < 2.5) and (abs(yd) < 2.5):
        numuPx = partNumuDetect.p()[1][0]
        numuPy = partNumuDetect.p()[1][1]
        numuPz = partNumuDetect.p()[1][2]
        numuE = math.sqrt(numuPx*numuPx + numuPy*numuPy + numuPz*numuPz )
        eNumuData.Fill(numuE)
    if (partNueDetect.weight() > 0.0):
      xs = partNueProd.x()
      ys = partNueProd.y()
      zs = partNueProd.z()
      xd = partNueDetect.x()
      yd = partNueDetect.y()
      zd = partNueDetect.z()
      dist = math.sqrt((xd-xs)*(xd-xs) + (yd-ys)*(yd-ys) +(zd-zs)*(zd-zs))
      nueDist.Fill(dist)
      if (abs(xd) < 2.5) and (abs(yd) < 2.5):
        nuePx = partNueDetect.p()[1][0]
        nuePy = partNueDetect.p()[1][1]
        nuePz = partNueDetect.p()[1][2]
        nueE = math.sqrt(nuePx*nuePx + nuePy*nuePy + nuePz*nuePz )
        eNueData.Fill(nueE)
# Flash neutrinos - if only transfer line decays of the pion are included
#    if (partNumuDetect.weight() > -1.0):
    if (partNumuDetect.weight() > 0.0):
      xs = partFlshNu.x()
      ys = partFlshNu.y()
      zs = partFlshNu.z()
      xd = partNumuDetect.x()
      yd = partNumuDetect.y()
      zd = partNumuDetect.z()
      dist = math.sqrt((xd-xs)*(xd-xs) + (yd-ys)*(yd-ys) +(zd-zs)*(zd-zs))
      flashNumuDist.Fill(dist)
      flashNumuSrcX.Fill(xs)
      flashNumuSrcY.Fill(ys)
      flashNumuSrcZ.Fill(zs)
      p = partFlshNu.p()
      flshNuSrcPX.Fill(p[1][0])
      flshNuSrcPY.Fill(p[1][1])
      flshNuSrcPZ.Fill(p[1][2])
      arrivalT.Fill(partNumuDetect.t())

histFile = os.path.join(StudyDir, ctrlInst.studyName())
print ("about to do ---- hm.histdo(histFile)")
hm.histdo(histFile)
print ("done ---- hm.histdo(histFile)")
fileName = os.path.join(StudyDir, ctrlInst.studyName() + "/plots" + str(ctrlInst.runNumber()) + ".root")
print (fileName)
print ("about to do ---- hm.histOutRoot(histFile)")
hm.histOutRoot(fileName)
fileName01 = os.path.join(StudyDir, ctrlInst.studyName() + "/nuDetector" + str(ctrlInst.runNumber()) + ".root")
hmDataOut.histOutRoot(fileName01)
print ("done ---- hm.histOutRoot(histFile)")

fileName02 = os.path.join(StudyDir, ctrlInst.studyName() + "/plots" + str(ctrlInst.runNumber()) + ".tex")
hm.texCreate(fileName02)
fileName03 = os.path.join(StudyDir, ctrlInst.studyName() + "/summary" + str(ctrlInst.runNumber()) + ".tex")
hC.summary(fileName03)

print()
print("========  analysing the eventHistory ends ========")
logging.info("========  analysing the eventHistory: end  ========")


sys.exit(0)
