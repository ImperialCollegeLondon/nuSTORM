#!/usr/bin/env python3

"""
Creates a shell script that can be run to submit all of
the Fluka simulations on FermiGrid computers,
using submitFermiGridJob.py for each job.

"""

import os
import sys

def run(inName, nJobs, nEvents, beamE, hornI, hornPol):

    runFileName = 'runJobs_{0}_{1:.0f}GeV_HPos.sh'.format(inName, beamE)
    if hornPol < 0:
        runFileName = runFileName.replace('HPos', 'HNeg')
    runFile = open(runFileName, 'w')
    
    for i in range(nJobs):
        line = 'python submitFermiGridJob.py {0} {1} {2} {3} {4} {5}\n'.format(inName, \
                                                     i+1, nEvents, beamE, hornI, hornPol)
        runFile.write(line)
        # Wait a few seconds
        runFile.write('sleep 2\n')
    
    runFile.close()
    
    
if __name__ == '__main__':

    inName = 'nuSTORMTarget'
    # 10k events take ~ 2 hours, 50k ~ 10 hours
    # 1 million events = 10k x 100 jobs or
    # 10 million events = 50k x 200 jobs (10 hrs each)
    nJobs = 100
    nEvents = 10000
    # Beam energy (GeV)
    beamE = 100.0
    # Horn current and polarity
    hornI = 230.0
    hornPol = 1.0
 
    nArg = len(os.sys.argv)
    if nArg > 1:
        inName = os.sys.argv[1]
    if nArg > 2:
        nJobs = int(os.sys.argv[2])
    if nArg > 3:
        nEvents = int(os.sys.argv[3])
    if nArg > 4:
        beamE = float(os.sys.argv[4])
    if nArg > 5:
        hornI = float(os.sys.argv[5])
    if nArg > 6:
        hornPol = float(os.sys.argv[6])
        
    run(inName, nJobs, nEvents, beamE, hornI, hornPol)
    
