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

nuSIMPATH = os.getenv('nuSIMPATH')
filename  = os.path.join(nuSIMPATH, \
                         '11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
nuStrt = nuPrdStrt.nuSTORMPrdStrght(filename)
                         

mc = muCnst.MuonConst()
pc = piCnst.PionConst()

nuEI = []
piEI = []

runNum = 1
eventNum = 1
#Not sure about the values


#Running PiEvt and NuEvt together
for i in range(10000):

    piEvt= piEvtInst.PionEventInstance(8)
    piEI.append(piEvt)
    tpi=piEvt.gettpi()
    piTraceSpaceCoord=piEvt.getTraceSpaceCoord()
    mu4mmtm=piEvt.getmu4mmtm()
    mucostheta=piEvt.getcostheta()
    nuEI.append(nuEvtInst.NeutrinoEventInstance(tpi,  piTraceSpaceCoord, mu4mmtm, mucostheta, nuStrt))

obj = eventHistory.eventHistory()
obj.outFile("testFile1.root")
obj.rootStructure()

#Storing the information in particle format
for (piEvt,nuEvt) in zip(piEI,nuEI):
  eventWeight = 1 
  #Target
  testParticle = particle.particle(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, "pi-")
  location = 'target'
  obj.addParticle(location, testParticle)
 
  #Production Straight 
  testParticle = particle.particle(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, "pi-")
  location = 'productionStraight'
  obj.addParticle(location, testParticle)
 
  # 1.Pion decay
  tpi = piEvt.gettpi()
  TSCpi = piEvt.getTraceSpaceCoord()
  ppiz = piEvt.getppiGen()
  testParticle = particle.particle(1, eventNum, TSCpi[0], TSCpi[1], TSCpi[2], TSCpi[3], TSCpi[4]*ppiz, TSCpi[5]*ppiz, ppiz, tpi, eventWeight, "pi-")
  location = 'pionDecay'
  obj.addParticle(location, testParticle)

  #Pion Flash
  numu4mmtm = piEvt.getnumu4mmtm()
  numu3mmtm = numu4mmtm[1]
  testParticle = particle.particle(1, eventNum, TSCpi[0], TSCpi[1], TSCpi[2], TSCpi[3], numu3mmtm[0], numu3mmtm[1], numu3mmtm[2], tpi, eventWeight, "numu")
  location="piFlashNu"
  obj.addParticle(location, testParticle)

   
  #Muon Production
  mu4mmtm = piEvt.getmu4mmtm()
  mu3mmtm = mu4mmtm[1]
  testParticle = particle.particle(1, eventNum, TSCpi[0], TSCpi[1], TSCpi[2], TSCpi[3], mu3mmtm[0], mu3mmtm[1], mu3mmtm[2], tpi, eventWeight, "mu-")
  location = "muonProduction"
  obj.addParticle(location, testParticle)

   
  if(nuEvt.getAbsorbed() == True):
     eventWeight=0
    
  #1. Muon Decay
  tmu = nuEvt.gettmu()
  TSCmu = nuEvt.getTraceSpaceCoord()
  pmuz = nuEvt.getpmu()
  testParticle = particle.particle(1, eventNum, TSCmu[0], TSCmu[1], TSCmu[2], TSCmu[3], TSCmu[4]*pmuz, TSCmu[5]*pmuz, pmuz, tmu, eventWeight, "mu-")
  location='muonDecay'
  obj.addParticle(location, testParticle)

   
  #2. Electron Production
  e4mmtm = nuEvt.gete4mmtm()
  e3mmtm = e4mmtm[1]
  testParticle = particle.particle(1, eventNum, TSCmu[0], TSCmu[1], TSCmu[2], TSCmu[3], e3mmtm[0], e3mmtm[1], e3mmtm[2], tmu, eventWeight, "e-")
  location='eProduction'
  obj.addParticle(location, testParticle)

  
  #3. nue Production
  nue4mmtm = nuEvt.getnue4mmtm()
  nue3mmtm = nue4mmtm[1]
  testParticle = particle.particle(1, eventNum, TSCmu[0], TSCmu[1], TSCmu[2], TSCmu[3], nue3mmtm[0], nue3mmtm[1], nue3mmtm[2], tmu, eventWeight, "nue")
  location='nueProduction'
  obj.addParticle(location, testParticle)
 
  
  #4. numu Production
  numu4mmtm = nuEvt.getnumu4mmtm()
  numu3mmtm = numu4mmtm[1]
  testParticle = particle.particle(runNum, eventNum, TSCmu[0], TSCmu[1], TSCmu[2], TSCmu[3], numu3mmtm[0], numu3mmtm[1], numu3mmtm[2], tmu, eventWeight, "numu")
  location='numuProduction'
  obj.addParticle(location, testParticle)
  
  # numu Detector
  testParticle = particle.particle(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, "numu")
  location = 'numuDetector'
  obj.addParticle(location, testParticle)
  
  # nue Detector
  testParticle = particle.particle(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, "nue")
  location = 'nueDetector'
  obj.addParticle(location, testParticle)
  obj.fill()
 
   
obj.write()
obj.outFileClose()   


