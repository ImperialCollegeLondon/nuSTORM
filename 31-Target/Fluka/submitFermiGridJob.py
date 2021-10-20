#!/usr/bin/env python3

"""
Python script to submit a Fluka simulation run on FermiGrid computers.

"""

import os
import platform
import random
import sys

class parameters(object):

    def __init__(self, inName, jobNum, nEvents, beamE, hornI, hornPol, useDPMJET = True):

        # Home directory
        self.home = os.getenv('HOME')
        # Fluka release and base directory, assumed to be in $HOME/app/FlukaArchive/flukaVer
        self.flukaVer = 'fluka4-1.1'
        self.flukaDir = '{0}/app/FlukaArchive/{1}'.format(self.home, self.flukaVer)
        # Location of flair and the working sub-dir containing input files
        self.flairVer = 'flair-3.1'
        self.flairDir = '{0}/{1}'.format(self.flukaDir, self.flairVer)
        self.workName = 'nuSTORM'
        # Location of input files
        self.inDir = '{0}/{1}'.format(self.flairDir, self.workName)

        # Horn B field parameters: current and polarity. The Fortran code is expecting
        # FieldPars.dat, so we use either the pos or negative file, as required
        self.fieldPars = 'FieldPars.dat'
        self.fieldParsPos = 'FieldParsPos.dat'
        self.fieldParsNeg = 'FieldParsNeg.dat'
        
        # Scratch area for batch jobs:
        self.user = os.getenv('USER')
        self.scratch = '/pnfs/dune/scratch/users/{0}/FlukaArchive'.format(self.user)
        self.baseJob = '{0}/{1}/{2}'.format(self.scratch, self.flukaVer, self.workName)
        # Name of Fluka "release" tar ball
        self.flukaTarName = 'FlukaInstall.tar.gz'
        self.flukaTarFile = '{0}/{1}'.format(self.baseJob, self.flukaTarName)

        # Proton beam energy (GeV)
        self.beamE = beamE
        # Beam FWHM size (cm)
        self.beamFWHM = 0.629
        
        # Horn current and polarity
        self.hornI = hornI
        self.hornPol = hornPol
        
        # Fluka executable
        self.flukaExe = 'flukahp'
        if useDPMJET == True:
            self.flukaExe = 'flukadpm3'
            
        # Fluka input template file
        self.inName = inName
        # Label including beam E and horn polarity
        self.label = '{0}_{1:.0f}GeV_HPos'.format(inName, beamE)
        if self.hornPol < 0.0:
            self.label = '{0}_{1:.0f}GeV_HNeg'.format(inName, beamE)
        
        # Job index number
        self.jobNum = jobNum
        # Number of events to simulate
        self.nEvents = nEvents
        
        # Name of the job and its scratch directory location
        self.jobName = '{0}_{1}'.format(self.label, self.jobNum)
        self.jobDir = '{0}/{1}'.format(self.baseJob, self.jobName)
        
        # Batch machine temporary directory location
        self.batchDir = '$_CONDOR_SCRATCH_DIR'
        
        # Operating system
        self.OSVer ="SL%i"%int(float(platform.linux_distribution()[1]))
        
        
def run(pars):

    print('Base job directory = {0}'.format(pars.baseJob))
    print('inName = {0}, label = {1}, job no = {2}, nEvents = {3}'.format(pars.inName, pars.label,
                                                                          pars.jobNum, pars.nEvents))

    # Run directory for the job
    print('JobName = {0}'.format(pars.jobName))
    print('JobDir = {0}'.format(pars.jobDir))
    if not os.path.exists(pars.jobDir):
        print('Creating {0}'.format(pars.jobDir))
        os.makedirs(pars.jobDir)

        os.chmod(pars.jobDir, 0777)
        
    # Create tar file containing Fluka release and copy it to scratch area
    if not os.path.exists(pars.flukaTarFile):
        print('Creating {0}'.format(pars.flukaTarFile))
        # Set the release as the current directory, leaving out flair as well as the
        # lower-level user directories (otherwise they would be needed upon extraction)
        tarCmd = 'tar --exclude={0} -czf {1} -C {2} .'.format(pars.flairVer, pars.flukaTarFile,
                                                              pars.flukaDir, pars.flukaTarName)
        print('tarCmd = {0}'.format(tarCmd))
        os.system(tarCmd)
    else:
        # Tar file has already been created
        print('Using tarball {0}'.format(pars.flukaTarFile))

    # Set horn A field parameter file and Fluka input file
    (inputFile, fieldFile) = createFlukaInput(pars)
        
    # Now create the script containing the job commands
    jobScript = createJobScript(pars, inputFile, fieldFile)

    # Finally submit the job; for tests, use --timeout=Xm to stop after X minutes (360s for N = 10)
    # 10k events take ~ 2 hours, 50k ~ 10 hours => 10 million events = 200 jobs. Use "12h"
    logFile = '{0}/submitJob.log'.format(pars.jobDir)
    jobCmd = 'jobsub_submit.py -N 1 --resource-provides=usage_model=OPPORTUNISTIC ' \
             '--expected-lifetime 12h --OS={0} --group=dune ' \
             '-L {1} file://{2}'.format(pars.OSVer, logFile, jobScript)
    print('jobCmd = {0}'.format(jobCmd))

    os.system(jobCmd)
    print('DONE\n')
    
    
def createFlukaInput(pars):

    # Copy Fluka input file, changing random seed given job number
    oldFileName = '{0}/{1}.inp'.format(pars.inDir, pars.inName)
    newFileName = '{0}/{1}.inp'.format(pars.jobDir, pars.label)
    print('Creating {0}\nbased on {1}'.format(newFileName, oldFileName))
    
    if os.path.exists(newFileName):
	os.remove(newFileName)
    newFile = open(newFileName, 'w')

    # Random number using 8 digit integers
    # Seed set using input file name and job number
    random.seed(a=oldFileName+str(pars.jobNum))
    randInt = random.randint(10000000, 99999999)
  
    with open(oldFileName, 'r') as f:
        for line in f:
            if 'RANDOMIZ' in line:
                newFile.write('RANDOMIZ   , {0:.1f}, {1:.1f}\n'.format(1, randInt))
            elif 'START' in line:
                # Number of events
                newFile.write('START      , {0:.1f}\n'.format(nEvents))
            elif 'BEAM' in line and 'BEAMPOS' not in line:
                # Beam energy and size
                line = '{0:<10}{1:>10.1f}{2:>10}{3:>10}{4:>10.3f}{5:>10.3f}{6:>16}\n'.format('BEAM', \
                                    -pars.beamE, 0.0, 0.0, -pars.beamFWHM, -pars.beamFWHM, 'PROTON')
                newFile.write(line)
            else:
                # Copy line unchanged
                newFile.write(line)
    
    newFile.close()

    # Copy horn field parameters to scratch job directory
    fieldFile = '{0}/{1}'.format(pars.inDir, pars.fieldParsPos)
    if pars.hornPol < 0:
        fieldFile = '{0}/{1}'.format(pars.inDir, pars.fieldParsNeg)

    fFile = open(fieldFile, 'w')
    fFile.write('{0:.1f} {1:.1f}\n'.format(pars.hornI, pars.hornPol))
    fFile.close()
    
    fieldCopy = '{0}/{1}'.format(pars.jobDir, pars.fieldPars)
    print('Copying {0}\nto {1}'.format(fieldFile, fieldCopy))
    os.system('cp -f {0} {1}'.format(fieldFile, fieldCopy))

    # Return copied filenames
    return (newFileName, fieldCopy)


def createJobScript(pars, inputFile, fieldFile):

    jobScript = '{0}/job{1}.sh'.format(pars.jobDir, pars.jobName)
    print('Creating job script {0}'.format(jobScript))
    if os.path.exists(jobScript):
	os.remove(jobScript)
    jobFile = open(jobScript, 'w')

    # Setup job environment and compiler
    jobFile.write('date\n')
    jobFile.write('source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups\n')
    jobFile.write('setup ifdhc\n')
    jobFile.write('source /cvmfs/sft.cern.ch/lcg/releases/gcc/9.2.0/x86_64-slc6/setup.sh\n\n')
    # Limit copying attempts to avoid stalled jobs
    jobFile.write('export IFDH_CP_MAXRETRIES=2\n')
    
    jobFile.write('echo Batch dir is {0}'.format(pars.batchDir))
    
    # Copy Fluka release tarball from job pnfs area to batch machine local scratch area
    batchTarFile = '{0}/{1}.tar.gz'.format(pars.batchDir, pars.jobName)
    jobFile.write('echo HereA: $_CONDOR_SCRATCH_DIR\n')
    jobFile.write('ifdh cp {0} {1}\n'.format(pars.flukaTarFile, batchTarFile))
    jobFile.write('echo HereB: $_CONDOR_SCRATCH_DIR\n')        
    jobFile.write('tar -xzf {0} -C {1}\n\n'.format(batchTarFile, pars.batchDir))

    jobFile.write('echo Here1: $_CONDOR_SCRATCH_DIR\n')
    
    # Setup Fluka environment
    jobFile.write('export FLUPRO={0}\n'.format(pars.batchDir))
    jobFile.write('export FLUFOR=gfortran\n')
    jobFile.write('export GFORFLU=gfortran\n')

    jobFile.write('echo FLUPRO = $_CONDOR_SCRATCH_DIR\n')
    jobFile.write('echo Here2: $_CONDOR_SCRATCH_DIR\n')

    # Copy Fluka input & field parameter files
    jobFile.write('ifdh cp {0} {1}/{2}\n'.format(inputFile, pars.batchDir, os.path.basename(inputFile)))
    jobFile.write('ifdh cp {0} {1}/{2}\n\n'.format(fieldFile, pars.batchDir, os.path.basename(fieldFile)))

    jobFile.write('echo Here3: $_CONDOR_SCRATCH_DIR\n')
    jobFile.write('ls -lrta {0}\n\n'.format(pars.batchDir))

    # Fluka run command
    jobFile.write('cd {0}\n'.format(pars.batchDir))
    runCmd = '$FLUPRO/bin/rfluka -e $FLUPRO/bin/{0} -N0 -M1 {1} \n\n'.format(pars.flukaExe,
                                                                             pars.label)
    jobFile.write(runCmd)

    jobFile.write('echo Here4\n')

    # List all output (to make sure everything was generated OK)
    jobFile.write('ls -lrt {0}\n\n'.format(pars.batchDir))
    
    # Copy output files. First create tarball of all usrbin histograms.
    # Scratch directory is usually /storage/local/data1/condor/execute/dir_X/no_xfer (X = integer)
    getBinList = 'ls {0}*fort* > binList.txt\n'.format(pars.label)
    jobFile.write(getBinList)
    # Then add event-by-event energy deposition and aperture tracking files
    jobFile.write('ls {0}*.txt >> binList.txt\n'.format(pars.label))
    # Finally add the log files
    jobFile.write('ls {0}*.log >> binList.txt\n'.format(pars.label))
    jobFile.write('ls {0}*.out >> binList.txt\n'.format(pars.label))
    jobFile.write('ls {0}*.err >> binList.txt\n'.format(pars.label))    
    binTarCmd = 'tar czf binList.tar.gz -T binList.txt\n'
    jobFile.write(binTarCmd)

    jobFile.write('echo Here5\n')

    # List all output (to make sure everything was generated OK)
    jobFile.write('ls -lrt {0}\n\n'.format(pars.batchDir))

    # Copy usrbin tarball. We can only extract it after the job has finished, since
    # dCache files are immutable within a job and can't be modified once copied
    jobFile.write('ifdh cp -D {0}/binList.tar.gz {1}/\n'.format(pars.batchDir, pars.jobDir))

    # Force the job to stop at the end
    jobFile.write('echo Copied all files\n')
    jobFile.write('exit 0\n')
    
    jobFile.close()

    # Make job script executable
    os.chmod(jobScript, 0777)

    return jobScript
    
    
if __name__ == '__main__':

    inName = 'nuSTORMTarget'
    jobNum = 1
    nEvents = 10
    beamE = 100.0 # or 26 GeV
    hornI = 230.0 # kA
    hornPol = 1.0 # +1 or -1
    
    nArg = len(os.sys.argv)
    if nArg > 1:
        inName = os.sys.argv[1]
    if nArg > 2:
        jobNum = int(os.sys.argv[2])
    if nArg > 3:
        nEvents = int(os.sys.argv[3])
    if nArg > 4:
        beamE = float(os.sys.argv[4])
    if nArg > 5:
        hornI = float(os.sys.argv[5])
    if nArg > 6:
        hornPol = float(os.sys.argv[6])
        
    useDPMJET = True
    pars = parameters(inName, jobNum, nEvents, beamE, hornI, hornPol, useDPMJET)
    run(pars)
