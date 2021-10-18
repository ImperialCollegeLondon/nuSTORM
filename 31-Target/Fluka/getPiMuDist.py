#!/usr/bin/env python3

"""
Python script to store the pion and muon tracking information,
from the Fluka BXDRAW subroutine output at the downstream end
of the horn, in ROOT format.

"""

import math
import os
import ROOT
import sys
from cppyy.ll import cast

ROOT.gInterpreter.Declare("""
struct FlukaStruct{
   int iEvt;
   int PDGId;
   double totE;
   double KE;
   double x;
   double y;
   double z;
   double rxy;
   double px;
   double py;
   double pz;
   double pT;
   double p;
   double tns;
   double weight;
};
""")

# Base directory location for the Fluka BXDRAW output files
baseDir = 'fluka4-1.1/nuSTORM'

def run(label, NJobs, NEvtJob):

    print('Label {0} for {1} jobs; {2} events per job'.format(label, NJobs, NEvtJob))

    rootFile = '{0}/{1}All/{1}.root'.format(baseDir, label)
    print('Creating {0}'.format(rootFile))
    
    outFile = ROOT.TFile.Open(rootFile, 'recreate')
    outFile.cd()
    data = ROOT.TTree('Data', 'Data')
    data.SetDirectory(outFile)

    fS = ROOT.FlukaStruct()
    
    data.Branch('iEvt', cast['void*'](ROOT.addressof(fS, 'iEvt')), 'iEvt/I')
    data.Branch('PDGId', cast['void*'](ROOT.addressof(fS, 'PDGId')), 'PDGId/I')
    data.Branch('totE', cast['void*'](ROOT.addressof(fS, 'totE')), 'totE/D')
    data.Branch('KE', cast['void*'](ROOT.addressof(fS, 'KE')), 'KE/D')
    data.Branch('x', cast['void*'](ROOT.addressof(fS, 'x')), 'x/D')
    data.Branch('y', cast['void*'](ROOT.addressof(fS, 'y')), 'y/D')
    data.Branch('z', cast['void*'](ROOT.addressof(fS, 'z')), 'z/D')
    data.Branch('rxy', cast['void*'](ROOT.addressof(fS, 'rxy')), 'rxy/D')
    data.Branch('px', cast['void*'](ROOT.addressof(fS, 'px')), 'px/D')
    data.Branch('py', cast['void*'](ROOT.addressof(fS, 'py')), 'py/D')
    data.Branch('pz', cast['void*'](ROOT.addressof(fS, 'pz')), 'pz/D')
    data.Branch('pT', cast['void*'](ROOT.addressof(fS, 'pT')), 'pT/D')
    data.Branch('p', cast['void*'](ROOT.addressof(fS, 'p')), 'p/D')
    data.Branch('tns', cast['void*'](ROOT.addressof(fS, 'tns')), 'tns/D')
    data.Branch('weight', cast['void*'](ROOT.addressof(fS, 'weight')), 'weight/D')

    # Each job should have the same number of events    
    for iJob in range(1, NJobs+1):

        dirName = '{0}/{1}_{2}'.format(baseDir, label, iJob)
        textFile = '{0}/{1}001_output.txt'.format(dirName, label)
        print('Processing {0}'.format(textFile))

        # Open tracking plane data text file
        with open(textFile, 'r') as inFile:
            
            for line in inFile:

                # Store space-separated line entries in array
                arr = line.split()

                if (len(arr) < 12):
                    print('Expecting 12 words\n')
                    continue

                # Event number and PDGId
                fS.iEvt = int(arr[0]) + (iJob-1)*NEvtJob
                partId = int(arr[1])
                fS.PDGId = getPDGId(partId)

                # Total energy and KE (GeV)
                fS.totE = float(arr[2])
                fS.KE = float(arr[3])

                # Weight and position (cm)
                fS.weight = float(arr[4])
                fS.x = float(arr[5])
                fS.y = float(arr[6])
                fS.z = float(arr[7])
                fS.rxy = math.sqrt(fS.x*fS.x + fS.y*fS.y)
                
                # Direction cosines
                cosx = float(arr[8])
                cosy = float(arr[9])
                cosz = float(arr[10])
                
                # Momentum (GeV)
                pSq = fS.KE*(2.0*fS.totE - fS.KE)
                if pSq > 0.0:
                    fS.p = math.sqrt(pSq)
                else:
                    fS.p = 0.0

                fS.px = fS.p*cosx
                fS.py = fS.p*cosy
                fS.pz = fS.p*cosz
                fS.pT = math.sqrt(fS.px*fS.px + fS.py*fS.py)
                
                # Lifetime (nanosecs)
                fS.tns = float(arr[11])*1e9

                # Store data
                data.Fill()

    # Write data and close output
    print('Writing {0}'.format(rootFile))
    data.Write()
    outFile.Close()

    
def getPDGId(partId):

    # Get the PDG integer from the Fluka id
    pdgId = 0
    
    if partId == 1:
        pdgId = 2212 # proton
    elif partId == 8:
        pdgId = 2112 # neutron
    elif partId == 13:
        pdgId = 211 # pi+
    elif partId == 14:
        pdgId = -211 # pi-
    elif partId == 23:
        pdgId = 111 # pi0
    elif partId == 15:
        pdgId = 321 # K+
    elif partId == 16:
        pdgId = -321 # K-
    elif partId == 19 or partId == 24 or partId == 25:
        pdgId = 310 # K0s, K0 or K0bar
    elif partId == 12:
        pdgId = 130 # K0L
    elif partId == 10:
        pdgId = -13 # mu+
    elif partId == 11:
        pdgId = 13 # mu-
    else:
        print('Using PDGId = {0} for PID partId = {1}'.format(pdgId, partId))
        
    return pdgId


def getMass(PDGId):

    # Get the particle rest mass for important PDGId's
    mass = 0.0

    if (PDGId == 2212):
        mass = 0.938272 # proton

    elif (PDGId == 2112):
        mass = 0.939565 # neutron

    elif (abs(PDGId) == 211):
        mass = 0.139570 # pi+-

    elif (PDGId == 111):
        mass = 0.134977 # pi0

    elif (abs(PDGId) == 321):
        mass = 0.493677 # K+-

    elif (PDGId == 310 or PDGId == 130):
        mass = 0.497614 # K0

    elif (abs(PDGId) == 13):
        mass = 0.105658 # mu+-

    else:
        print('getMass: Ignored PDGId {0}'.format(PDGId))

    return mass
  

if __name__ == '__main__':

    label = 'nuSTORMTarget_26GeV_HPos'
    NJobs = 100
    NEvtJob = 10000 # Number of events per job
    
    NArg = len(sys.argv)
    
    if NArg > 1:
        label = sys.argv[1]
    if NArg > 2:
        NJobs = int(sys.argv[2])
    if NArg > 3:
        NEvtJob = int(sys.argv[3])
        
    run(label, NJobs, NEvtJob)
