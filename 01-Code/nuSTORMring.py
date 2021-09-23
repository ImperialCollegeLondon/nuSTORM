import os
import MuonConst as muCnst
import PionConst as piCnst
import nuSTORMPrdStrght as nuPrdStrt
import NeutrinoEventInstance as nuEvtInst
import PionEventInstance as piEvtInst
import numpy as np
import Simulation as Simu
import sys
import particle as particle
import eventHistory as eventHistory

mc = muCnst.MuonConst()
pc = piCnst.PionConst()

nuEI = []
piEI = []
obj = eventHistory.eventHistory()

runNum = 0
eventNum = 0
eventWeight = 0 #Not sure about the values
pdgCode = 0
electronmass = 9.11*10**-31 # not sure about the units

#Running PiEvt and NuEvt together
for i in range(1000):

    piEvt= piEvtInst.PionEventInstance(8)
    piEI.append(piEvt)
    tpi=piEvt.gettpi()
    piTraceSpaceCoord=piEvt.getTraceSpaceCoord()
    mu4mmtm=piEvt.getmu4mmtm()
    nuEI.append(nuEvtInst.NeutrinoEventInstance(tpi,  piTraceSpaceCoord, mu4mmtm, filename=1))

#Deleting the muons that are not accepted:
for nuEvt in nuEI:
    if(nuEvt.getAbsorbed() == True):
      nuEI.remove(nuEvt) 


#Storing the information in particle format
for piEvt in piEI:
   # 1.Pion decay
   tpi = piEvt.gettpi()
   TSCpi = piEvt.getTraceSpaceCoord()
   ppiz = piEvt.getppiGen()
   testParticle = particle.particle(runNum, eventNum, TSCpi[1], TSCpi[2], TSCpi[3], TSCpi[0], TSCpi[4]*ppiz, TSCpi[5]*ppiz, ppiz, tpi, eventWeight, pc.mass(), pdgCode)
   location = 'pionDecay'
   obj.addParticle(location, testParticle)
   retPar = obj.findParticle(location)
   
   #Pion Flash
   numu4mmtm = piEvt.getnumu4mmtm()
   numu3mmtm = numu4mmtm[1]
   testParticle = particle.particle(runNum, eventNum, TSCpi[1], TSCpi[2], TSCpi[3], TSCpi[0], numu3mmtm[0], numu3mmtm[1], numu3mmtm[2], tpi, eventWeight, 0, pdgCode)
   location="piFlashNu"
   obj.addParticle(location, testParticle)
   retPar = obj.findParticle(location)
   
   #Muon Production
   mu4mmtm = piEvt.getmu4mmtm()
   mu3mmtm = mu4mmtm[1]
   testParticle = particle.particle(runNum, eventNum, TSCpi[1], TSCpi[2], TSCpi[3], TSCpi[0], mu3mmtm[0], mu3mmtm[1], mu3mmtm[2], tpi, eventWeight, mc.mass(), pdgCode)
   location = "muonProduction"
   obj.addParticle(location, testParticle)
   retPar = obj.findParticle(location)
   
for nuEvt in nuEI:
  #1. Muon Decay
  tmu = nuEvt.gettmu()
  TSCmu = nuEvt.getTraceSpaceCoord()
  pmuz = nuEvt.getpmu()
  testParticle = particle.particle(runNum, eventNum, TSCmu[1], TSCmu[2], TSCmu[3], TSCmu[0], TSCmu[4]*pmuz, TSCmu[5]*pmuz, pmuz, tmu, eventWeight, mc.mass(), pdgCode)
  location='muonDecay'
  obj.addParticle(location, testParticle)
  retPar = obj.findParticle(location)
   
  #2. Electron Production
  e4mmtm = nuEvt.gete4mmtm()
  e3mmtm = e4mmtm[1]
  testParticle = particle.particle(runNum, eventNum, TSCmu[1], TSCmu[2], TSCmu[3], TSCmu[0], e3mmtm[0], e3mmtm[1], e3mmtm[2], tmu, eventWeight, electronmass, pdgCode)
  location='eProduction'
  obj.addParticle(location, testParticle)
  retPar = obj.findParticle(location)
  
  #3. nue Production
  nue4mmtm = nuEvt.getnue4mmtm()
  nue3mmtm = nue4mmtm[1]
  testParticle = particle.particle(runNum, eventNum, TSCpi[1], TSCpi[2], TSCpi[3], TSCpi[0], nue3mmtm[0], nue3mmtm[1], nue3mmtm[2], tpi, eventWeight, 0, pdgCode)
  location='nueProduction'
  obj.addParticle(location, testParticle)
  retPar = obj.findParticle(location)
  
  #4. numu Production
  numu4mmtm = nuEvt.getnumu4mmtm()
  numu3mmtm = numu4mmtm[1]
  testParticle = particle.particle(runNum, eventNum, TSCpi[1], TSCpi[2], TSCpi[3], TSCpi[0], numu3mmtm[0], numu3mmtm[1], numu3mmtm[2], tpi, eventWeight, 0, pdgCode)
  location='numuProduction'
  obj.addParticle(location, testParticle)
  retPar = obj.findParticle(location)
    
obj.display()
