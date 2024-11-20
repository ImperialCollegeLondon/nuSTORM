#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for going through the eventHistory and checking simulation
=================================================================

  Assumes that nuSim code is in python path.

    @file checkGeneration.py

    @brief   anaylyse the eventHistory

    @author  Paul Kyberd

    @version     1.5
    @date        27 October 2024

    Add piFlashNu to stations

    @version     1.4
    @date        26 October 2024

    Add numuDetector to stations

    @version     1.3
    @date        25 October 2024

    Add a parameter to set a particle

    @version     1.2
    @date        25 October 2024

    Add a weight cut parameter to allow plots with or without weight cut
    Set to true by default

    @version     1.1
    @date        25 October 2024

    add an explanation of its operation
    make the station selectable by passing a parameter of answering a prompt

    @version     1.0
    @date        23 October 2024


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
import historyChecks as hChecks
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
helpStr='Set a number of key-VALUE PAIRS:some stuff\n'
helpStr= helpStr + 'some more stuff\n'
parser.add_argument("--set", metavar="KEY=VALUE", 
                    nargs='+', 
                    help='Set a number of key-VALUE PAIRS: --set explain=True for help', 
                    action=keyvalue)
args = parser.parse_args()
print(args.set)
#  set flags only if the parameters dictionary was created - otherwise no parameters
runNumberFlg = False
cntrlFlg = False
dbgSetFlg = False
statiionFlg = False
wtCutSetFlg = False
particleSetFlg = False
if args.set:
  if "explain" in args.set.keys():
    print ("Recognised keys are: explain; runNumber; cntrl; dbg; station")
    print ("explain=True                this text")
    print ("station=stationName         to check a particular station ")
    print ("   Allowed values are:      Target; PSStart; PSEnd; piFlashNu; muonDecay")
    print ("runNumber=<number>")
    print ("wtCut=False                 Cut on the weight. Must be greater than 1 - normally on, switches off" )
    print ("dbg=True                    False also allowed, but is the default")
    print ("particle=<part>             part can be muon or pion")
    exit()
  if "runNumber" in args.set.keys():
    print ("runNumber  in dict")
    runNumberFlg = True 
    runNumberVal = args.set["runNumber"]
  else:
    print ("runNumber not in dict")

  if "cntrl" in args.set.keys():
    print ("run control file in args dict")
    cntrlFlg = True 
    cntrlFile = args.set["cntrl"]
  else:
    print ("cntrl file not in args.dict")

  if "station" in args.set.keys():
    print ("station in args dict")
    stationFlg = True 
    stationName = args.set["station"]
  else:
    print ("cntrl file not in args.dict")

  if "dbg" in args.set.keys():
    print ("dbg in args dict")
    dbgSetFlg = True 
    dbgFlg = args.set["dbg"]
  else:
    print ("dbg not in args.dict")

  if "wtCut" in args.set.keys():
    print ("wtCut in args dict")
    wtCutSetFlg = True 
    wtCutFlg = args.set["wtCut"]
  else:
    print ("wtCut not in args.dict")
    wtCutFlg = True


  if "particle" in args.set.keys():
    print ("particle in args dict")
    particleSetFlg = True 
    particleType = args.set["particle"]
    if particleType == "muon": 
      particleID = -13
    elif particleType == "pion": 
      particleID = 211
    else:
      print ("particle type is ", particleID, "  not recognised - muon or pion allowed")
  else:
    particleID = 0
    particleType = "all"


#   Set debug flag - if nothing default to false
  if not dbgSetFlg:
    dbgFlg = False

#   Set weight cut flag - if nothing default to false
#  print (" 1 --- wtCutSetFlg is ", wtCutSetFlg, "   and wtCutFlg is ", wtCutFlg)
#  if not wtCutSetFlg:
#    wtCutFlg = True
#  print (" 2 --- wtCutSetFlg is ", wtCutSetFlg, "   and wtCutFlg is ", wtCutFlg)

#   Get environment variables
StudyDir = os.getenv('StudyDir')
StudyName = os.getenv('StudyName')
nuSIMPATH = os.getenv('nuSIMPATH')

#   Get name of control file
if not cntrlFlg:
  cntrlFile = input ("name of control file (exclude .dict) ")
controlFile = os.path.join(StudyDir, StudyName, cntrlFile + ".dict")
#   and open it
ctrlInst = control.control(controlFile)

#   Get name of station if not included
if not stationFlg:
  stationName = input ("name of station: Target; PSStart; PSEnd; numuDet;")
#   and set it 
if   stationName == "Target":       station = "Target"
elif stationName == "PSStart":      station = "productionStraightStart"
elif stationName == "PSEnd":        station = "productionStraightEnd"
elif stationName == "piFlashNu":    station = "piFlashNu"
elif stationName == "muonDecay":    station = "muonDecay"
elif stationName == "numuDetector": station = "numuDetector"
else:
  exit("checkStation - station name unrecognised " + stationName)

#   Open the input file - either the current run number or one passed in the parameters
#     and an associated output file/12-examp
if runNumberFlg:
  rootFilename = os.path.join(StudyDir, StudyName, 'normalisation' + str(runNumberVal) +'.root')
  rootOutFilename = os.path.join(StudyDir, StudyName, 'checkPlots' + str(runNumberVal) +'.root')
else:
  rootFilename = os.path.join(StudyDir, StudyName, 'normalisation' + str(ctrlInst.runNumber()) +'.root')
  rootOutFilename = os.path.join(StudyDir, StudyName, 'checkPlots' + str(ctrlInst.runNumber()) +'.root')
print(f"rootFileName is {rootFilename}")

#       logfile initialisation
logging.basicConfig(filename=ctrlInst.logFile(), encoding='utf-8', level=logging.INFO)
logging.info("\n\n")
logging.info("========  checks on the eventHistory: start  ======== Version %s", aHVersion)
logging.info("  Control file %s", controlFile)
logging.info("  root input file %s ", rootFilename)
logging.info("  station is %s", station)
logging.info("  weight cut is %s", wtCutFlg)
logging.info("  particleType is %s", particleType)

print("========  checks on the eventHistory start  ========")
print("  control file is ", controlFile)
print("  root input file is ", rootFilename)
print("  station is ", station)
print("  weight cut is ", wtCutFlg)
print("  particleType is ", particleType)

#   Create an eventHistory object which can be filled
objRd = eventHistory.eventHistory()
#   Connect it to the root file
objRd.inFile(rootFilename)

#   Get the number of events - and output them
nEvent = objRd.getEntries()
print ("  Number of entries is ", nEvent)
logging.info("  Number of events is  %s ", nEvent)

#   Set up the histogram manager
hm = histoManager.histoManager()

#   Get the dictionary containing the cuts
cutsDict = ctrlInst.plotsDict()

#   Set up the history checks object
hChck = hChecks.historyChecks(hm)
if station == "Target":
  hChck.createTargetPlots(wtCutFlg, cutsDict, particleID)
if station == "productionStraightStart":
  hChck.createPStStartPlots(wtCutFlg, cutsDict, particleID)
if station == "productionStraightEnd":
  hChck.createPSEndPlots(wtCutFlg, cutsDict, particleID)
if station == "muonDecay":
  hChck.createmuonDecayPlots(wtCutFlg, cutsDict, particleID)
if station == "piFlashNu":
  hChck.createPiFlashNuPlots(wtCutFlg, cutsDict, particleID)
if station == "numuDetector":
  hChck.createnumuDetectorPlots(wtCutFlg, cutsDict, particleID)


print("start reading")
for pnt in range(nEvent):
# read an event
    if (pnt%1000 == 0): print ("Event ", pnt)
    objRd.readNext()
    if dbgFlg: objRd.display()
    if station == "Target":
        partTar = objRd.findParticle("target")
        hChck.fillTargetPlots(partTar)
    if station == "productionStraightStart":
        partPSStart = objRd.findParticle("productionStraight")
        hChck.fillPSStartPlots(partPSStart)
    if station == "productionStraightEnd":
        partPSEnd = objRd.findParticle(("prodStraightEnd"))
        hChck.fillPSEndPlots(partPSEnd)
    if station == "muonDecay":
        partmuDcy = objRd.findParticle(("muonDecay"))
        hChck.fillmuonDecayPlots(partmuDcy)
    if station == "piFlashNu":
        partPiFlashNu = objRd.findParticle(("piFlashNu"))
        hChck.fillPiFlashNuPlots(partPiFlashNu)
    if station == "numuDetector":
        partnumuDet = objRd.findParticle(("numuDetector"))
        hChck.fillnumuDetectorPlots(partnumuDet)

hm.histOutRoot(rootOutFilename)

print("========  checking the eventHistory ends ========")
logging.info("========  checking the eventHistory: end  ========")


sys.exit(0)

