'''
    genPion.py

  Assumes that nuSim code is in python path.

    @file genPion.py

    @brief   Generate a pion according to various distribtuions

    @author  Paul Kyberd

    @version     1.0
    @date        07 October 2024

    Generate a pion according to a number of different distributions.

    The code returns the px, py, pz, x, y, xp, yp of the pion
        xp = px/pz
        yp = py/pz

    @version     1.1
    @date        22 October 2024

    Add logging to the code


    @version     1.2
    @date        02 December 2024

    Instantiation includes the name of the fluka file to write to


    Version 1.2                                02 December 2024
    Author:                                     Paul Kyberd

    @todo:          Sort out the filename for the root histograms so it is consistent
                    with the mean pion momentum


'''

import numpy as np
import scipy as sp
import traceSpace as tSp
import RandomGenerator as Rndm
import KDEdataFromFluka as KDEFluka
import logging
import ROOT

class genPion:

    _dbg = True

    def __init__(self, pionMeanP, genType, nEvents, inputFile, KDEFlukaFile):

        '''
        @parameter      pionMeanP       mean momentum of the pion - but the target inpout file is hardwired
        @parameter      genType         can be "pencil", "root", "fluka"
        @parameter      nEvents         a bit clunky only actually used for KDE, but keep it simple
        @parameter      inputFile       this is the file for the root or KDE generation as appropriate
        '''

        self._Version = 1.0
        self._px = 0.0
        self._py = 0.0
        self._pz = 0.0
        self._xp = 0.0
        self._yp = 0.0
        self._x  = 0.0
        self._y  = 0.0
        self._pionMeanP = pionMeanP
        self._nEvents = nEvents

#   set up logging
        logging.getLogger("nuSIM")
        logging.info("genPion - Pion generation: %s", genType)

        if genType == "pencil":
            self._genType = genType
        elif genType == "parabolic":
            self._genType = genType
        elif genType == "root":
            self._genType = genType
            geV = int(pionMeanP)                        # Integer part of mean momentum
#            fracE = int((pionMeanP*10 - geV*10))          # 2 sig fig
            fracE = int((pionMeanP*100 - geV*100))        # 3 sig fig
            if self._dbg: print (f"geV is {geV} and fracE is {fracE}")
            histName = 'histP'+str(geV)+str(fracE)+'GeV'
            histName2Dx = 'histXPS'+str(geV)+str(fracE)+'GeV'
            histName2Dy = 'histYPS'+str(geV)+str(fracE)+'GeV'
            rootInputFilename = inputFile
            if self._dbg: 
                print (f"histName is {histName}")
                print (f"histName2Dx is {histName2Dx} - histName2Dy is {histName2Dy}")
            self._RndmGen = Rndm.RandomGenerator(rootInputFilename,histName,histName2Dx,histName2Dy)
        elif genType == "flukaKDE":
            self._genType = genType
            nSamples = self._nEvents
        #   Instantiate the KDE class
            self._KDEdata = KDEFluka.KDEdataFromFluka()
        #   Read the fluka file
            flukaDataFile = inputFile
            self._df = self._KDEdata.readFile(flukaDataFile)
        #   Generate the pion kinematics and put in a pandas array. num_samples is used set to the new sample size (hardwired)
            self._KDEdist = self._KDEdata.beamGen_kde(self._df, columns=['x', 'y', 'xp', 'yp', 'pz', 'p'], num_samples=nSamples)
        #   Save the new data to a text file
            self._KDEFile = KDEFlukaFile
            self._KDEdata.writeFile(self._KDEdist,self._KDEFile)
        #   Read the generated data file
            self._KDEdata.openGenFile(self._KDEFile)
        else:
            exit("unrecognised generation type")

    def __repr__(self):
        return "genPion"

    def __str__(self):
        retString =  "genPion: pion Generator - Version " + str(self._Version) + "\n"\
        + "Generation type is " + str(self._genType) + "\n x, y, px, py, pz, xp, yp -- "\
        + str(self._x) + ", " + str(self._y) + ", " + str(self._px) + ", " + str(self._py) + ", "\
        + str(self._xp) + ", " + str(self._y)
    
        return retString
    
    def newParticle(self):

        if self._genType == "pencil":
            pionTSp, pPion  = self.pencilGen()
        elif self._genType == "parabolic":
            pionTSp, pPion = self.parabolicGen()
        elif self._genType == "root":
            pionTSp, pPion = self.rootGen()
        elif self._genType == "flukaKDE":
            pionTSp, pPion = self.flukaKDEGen()
        else:
            exit("generation type not implemented")
        
#        if (self._dbg): print ("newParticle: pionTSp ", pionTSp)

        return pionTSp, pPion

    def pencilGen(self):
# generate with s = 0.0, x = 0.0, y = 0.0, z = 0.0, xp =  0.0, yp = 0.0
        s = 0.0
        x = 0.0
        y = 0.0
        z = 0.0
        xp = 0.0
        yp = 0.0
        pion = tSp.traceSpace(s, x, y, z, xp, yp)
        return pion, self._pionMeanP

    def parabolicGen(self):
# generate with s = 0.0, x = parabolic, y = parabolic, z = 0.0, xp =  parabolic, yp = parabolic
        s = 0.0
        x = 0.0
        y = 0.0
        z = 0.0
        xp = 0.0
        yp = 0.0
        pion = tSp.traceSpace(s, x, y, z, xp, yp)
        return pion, self._pionMeanP            # this line is wrong - need the proper pion momentum not mean

    def rootGen(self):
# generate from the root distributions derived from the fluka files
        pPion = self._RndmGen.getRandom()
        x, xp = self._RndmGen.getRandom2Dx()
        y, yp = self._RndmGen.getRandom2Dy()
        pz = np.sqrt(pPion**2/(1+xp**2+yp**2))
        s = 0.0
        z = -50.0
        pion = tSp.traceSpace(s, x, y, z, xp, yp)
        return pion, pPion

    def flukaKDEGen(self):
        ''' generate from the KDE derived from the fluka fluxes '''
        pion, pPion = self._KDEdata.readGenFile()
        return pion, pPion

if __name__ == "__main__" :

    print("does nothing")
