'''
    genPionTst.py

  Assumes that nuSim code is in python path.

    @file genPionTst.py

    @brief   Test genPion.py

    @author  Paul Kyberd
    @version     1.0
    @date        07 October 2024

    @author  Paul Kyberd
    @version     1.1
    @date        15 October 2024

    Adding tests when generating with the fluka KDE.
    Problem trying to create histos for more than one running condition
    causes a problem. At present only run 1 set
   
    Version 1.0                                 07 October 2024
    Author:                                     Paul Kyberd
'''

import scipy as sp
import genPion as genPion
import traceSpace as tSp
import histoManager as histoManager
import RandomGenerator as Rndm
import sys

class genPionTst:

    def __init__(self):

        self._hm = histoManager.histoManager()


    def main(self):

##! Start:
        classToTest="pionGen"
        print()
        print("========  ", classToTest, ": tests start  ========")
        testFails = 0
        pionTest = 1

##! Create instance of a pencil beam and test built-in methods:
        print()
        print("pionGen:", pionTest, " Create pionGen('pencil') and print quantities.")    

        meanP = 3.0
        nEvents = 1000
        inputFile = ""
        gp = genPion.genPion(meanP, "pencil", nEvents, inputFile)
        print ("instance of genPion ", gp)

        piTSp,p = gp.newParticle()

        print ("pion tracespace is ", piTSp)

##! Generate events with a pencil beam
        print()
        pionTest = pionTest + 1
        print("pionGen:", pionTest, " Generate events.")    
##! book the histograms
#        self.histBook()
        meanP = 3.0
        gp = genPion.genPion(meanP, "pencil", nEvents,  inputFile)
#! generate events and plot them        
#        for eventNum in range(5000):
#            pionStart, p = gp.newParticle()
#            self.histFill(pionStart, p)

##! output the histograms
#        self.histWrite('genPionPencilTst.root')

##  Create an instance of the parabolic beam and generate an events
        pionTest = pionTest + 1
        print()

###################################################################################################
##! Generate an event with a pencil beam
        print("pionGen:", pionTest, " Create pionGen('parabolic') and print quantities.")    
        meanP = 5.0
        inputFile = ""
        gp = genPion.genPion(meanP, "parabolic", nEvents, inputFile)
        print ("instance of genPion ", gp)

        print()
        print("========  ", classToTest, ": tests complete  ========")

        print ("\nNumber of tests is ", pionTest, " number of fails is ", testFails)

###################################################################################################
##! Generate events according to the root histograms of the fluka data
        print("pionGen:", pionTest, " Create pionGen('root') and print quantities.")    
##! book the histograms
        meanP = 5.0
#        self.histBook(meanP)
#        gp = genPion.genPion(meanP, "root")
#! generate events and plot them        
#        for eventNum in range(5000):
#            pionStart, p = gp.newParticle()
#            self.histFill(pionStart, p)

##! output the histograms
#        self.histWrite('genPionRootTst.root')

###################################################################################################
##! Generate events according to the root histograms of the fluka data
        print("pionGen:", pionTest, " Create pionGen('fluka') and print quantities.")    
##! book the histograms
        meanP = 5.0
        nEvents = 100000
        self.histBook(meanP)
        flukaDataFile = '/Users/paulkyberd/work/nuSTORMAug24/01-nuSIM/31-Target/v2/selected_pions_5GeV_2.0mm.txt'
        gp = genPion.genPion(meanP, "flukaKDE", nEvents, flukaDataFile)
#! generate events and plot them        
        for eventNum in range(nEvents):
            pionStart, p = gp.newParticle()
            print (" pionStart is ", pionStart)
            print (" p ", p)
            self.histFill(pionStart, p)

##! output the histograms
        self.histWrite('genPionKDETst.root')

        sys.exit(0)


###################################################################################################
    def histBook(self, meanP):
##! book the histograms
        self._hm = histoManager.histoManager()
        hBins = 100

        hTitle = "x"
        hLower = 0.0
        hUpper = 0.0
        self.xPlt = self._hm.book(hTitle, hBins, hLower, hUpper)
        hTitle = "y"
        self.yPlt = self._hm.book(hTitle, hBins, hLower, hUpper)
        hTitle = "z"
        hLower = 0.0
        hUpper = 0.0
        self.zPlt = self._hm.book(hTitle, hBins, hLower, hUpper)
        hLower = 0.0
        hUpper = 0.0
        hTitle = "s"
        self.sPlt = self._hm.book(hTitle, hBins, hLower, hUpper)
        hLower = 0.0
        hUpper = 0.0
        hTitle = "Xp"
        self.XpPlt = self._hm.book(hTitle, hBins, hLower, hUpper)
        hTitle = "Yp"
        self.YpPlt = self._hm.book(hTitle, hBins, hLower, hUpper)
        hTitle = "p"
        hLower = 0.0
        hUpper = 0.0     
        self.pPlt = self._hm.book(hTitle, hBins, hLower, hUpper)
        hLr = 0.0
        hUr = 0.0  
        hLp = 0.0
        hUp = 0.0  
        hTitle = "x v y"
        self.xyPlt = self._hm.book2(hTitle, hBins, hLp, hUp, hBins, hLp, hUp)
        hTitle = "x v xp"
        self.xxpPlt = self._hm.book2(hTitle, hBins, hLp, hUp, hBins, hLr, hUr)
        hTitle = "y v yp"
        self.yypPlt = self._hm.book2(hTitle, hBins, hLp, hUp, hBins, hLr, hUr)
        hTitle = "xp v yp"
        self.xpypPlt = self._hm.book2(hTitle, hBins, hLr, hUr, hBins, hLr, hUr)
        
        return

    def histFill(self, pionStart, p):
##! fill the histograms
        self.sPlt.Fill(pionStart.s())
        self.xPlt.Fill(pionStart.x())
        self.yPlt.Fill(pionStart.y())
        self.zPlt.Fill(pionStart.z())
        self.XpPlt.Fill(pionStart.xp())
        self.YpPlt.Fill(pionStart.yp())
        self.pPlt.Fill(p)

        self.xyPlt.Fill(pionStart.x(), pionStart.y())
        self.xxpPlt.Fill(pionStart.x(), pionStart.xp())
        self.yypPlt.Fill(pionStart.y(), pionStart.yp())
        self.xpypPlt.Fill(pionStart.xp(), pionStart.yp())


    def histWrite(self, fileName):
##! output the histograms      
        self._hm.histOutRoot(fileName)

if __name__ == "__main__" :
    gpt = genPionTst()
    gpt.main()

