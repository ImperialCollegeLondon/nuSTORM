"""
Model for calculating normalised numbers
========================================

@date: 05 Jan 2022
@author: Paul Kyberd
@version: 2.0

@description: read in a the dictionary with the histogram values in from an external file

"""

from pathlib import Path
import json

class histsCreate():


    def __init__(self, hM):
        self._locsStart = []
        self._locs=[]
        self._hists=[]
        self._hm = hM
        self._nHists = 10
        self._count = 0
        self._targetWeight = 0
        self._targetNoWeight = 0
        self._NZWeight = {"target":0, "productionStraight": 0, "pionDecay": 0, "muonProduction": 0,"piFlashNu": 0,
                "muonDecay": 0, "eProduction": 0, "numuProduction":0, "nueProduction": 0, "numuDetector": 0, "nueDetector": 0}
        self._zeroWeight = {"target":0, "productionStraight": 0, "pionDecay": 0, "muonProduction": 0,"piFlashNu": 0,
                "muonDecay": 0, "eProduction": 0, "numuProduction":0, "nueProduction": 0, "numuDetector": 0, "nueDetector": 0}

#  read in the dictionaries
        plotsFile = "04-Studies/plotsOne.dict"
        with open(plotsFile) as plotsFile:
            self._plotsInfo = json.load(plotsFile)

    def __repr__(self):
        return "create a set of histograms for a particle at a location in the event History"

    def histAdd(self, eventType):
        self._locsStart.append(len(self._hists))
        self._locs.append(eventType)
#        print ("locStart is ", self._locsStart)
#        print ("locs ", self._locs)
        hTitle = eventType + ":x"
        hBins  = 100

        hLower = self._plotsInfo['xLower'][eventType]
        hUpper = self._plotsInfo['xHigher'][eventType]
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))

        hLower = self._plotsInfo['yLower'][eventType]
        hUpper = self._plotsInfo['yHigher'][eventType]
        hTitle = eventType + ":y"
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))

        hLower = self._plotsInfo['zLower'][eventType]
        hUpper = self._plotsInfo['zHigher'][eventType]
        hTitle = eventType + ":z"
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))


        hTitle = eventType + ":weight"
        hLower = -10.0
        hUpper = 100.0
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))

        hLower = self._plotsInfo["sLower"][eventType]
        hUpper = self._plotsInfo["sHigher"][eventType]
        hTitle = eventType + ":s"
        if eventType == "muonDecay" or eventType == "eProduction":
            hBins = 50
        else:
            hBins = 100
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))

        hLower = self._plotsInfo["tLower"][eventType]
        hUpper = self._plotsInfo["tHigher"][eventType]
        hTitle = eventType + ":t"
        if eventType == "muonDecay" or eventType == "eProduction":
            hBins = 50
        elif eventType == "target":
            hbins = 70
            hUpper = 70.
        else:
            hBins = 100
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))

        hBins = 100
        hTitle = eventType + ":px"
        hLower = self._plotsInfo["pxLower"][eventType]
        hUpper = self._plotsInfo["pxHigher"][eventType]
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))

        hTitle = eventType + ":py"
        hLower = self._plotsInfo["pyLower"][eventType]
        hUpper = self._plotsInfo["pyHigher"][eventType]
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))

        hTitle = eventType + ":pz"
        hLower = self._plotsInfo["pzLower"][eventType]
        hUpper = self._plotsInfo["pzHigher"][eventType]
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))

        hTitle = eventType + ":E_v_t"
        self._hists.append(self._hm.book2(hTitle, 50, 0.0, 5.0, 50, 10000.0, 11000.0))


#        print (self._hists)
        return

    def histsFill(self, location, particle):

#  Get the location
        pnt = self._locs.index(location)
        hPnt = self._locsStart[pnt]

        x = particle.x()
        y = particle.y()
        z = particle.z()
        px = particle.p()[1][0]
        py = particle.p()[1][1]
        pz = particle.p()[1][2]
        wt = particle.weight()
        s = particle.s()
        t = particle.t()

        if (wt !=0):
            n = self._NZWeight.get(location)
            self._NZWeight[location] = n + 1
        else:
            n = self._zeroWeight.get(location)
            self._zeroWeight[location] = n + 1

        if (wt > 0.0):
            self._hists[hPnt+3].Fill(wt)
            self._hists[hPnt+4].Fill(s)
            self._hists[hPnt+5].Fill(t)
            self._hists[hPnt+6].Fill(px)
            self._hists[hPnt+7].Fill(py)
            self._hists[hPnt+9].Fill(pz,t)

            self._hists[hPnt+0].Fill(x)
            self._hists[hPnt+1].Fill(y)
            self._hists[hPnt+2].Fill(z)

            if (location == 'productionStraight'):
# plots with non-zero weight - for production straight
                 self._count = self._count + 1

            if (location == 'numuDetector'):
                if ((abs(x) < 2.5) and (abs(y) < 2.5)):
                        self._hists[hPnt+8].Fill(pz)
            else:
                self._hists[hPnt+8].Fill(pz)



    def summary(self):

        texHline = "\\hline\n"
        slash = "\\"
        print(" Summary of results")
        print (self._zeroWeight)
        print (self._NZWeight)
#  Outout the information as a latex table
        sumDir = Path.cwd()
        sumFile = str(Path.cwd())+"/summary.tex"
        print (sumFile)
        f = open(sumFile, "w")
        f.write("\\begin{tabular}{|  l | l | l  | p{9.0cm} | }\n" )
        f.write(texHline )
        f.write("Location   &  zero weight  &  Full weight & comments \\\\ \n"  )
        f.write(texHline )
        f.write("target                   & " + str(self._zeroWeight['target']) + " & " + str(self._NZWeight['target']) + " &\\\\ \n")
        f.write("productionStraight       & " + str(self._zeroWeight['productionStraight']) + " & " + str(self._NZWeight['productionStraight']) + " &\\\\ \n")
        f.write("pionDecay                & " + str(self._zeroWeight['pionDecay']) + " & " + str(self._NZWeight['pionDecay']) + " &\\\\ \n")
        f.write("muonProduction           & " + str(self._zeroWeight['muonProduction']) + " & " + str(self._NZWeight['muonProduction']) + " &\\\\ \n")
        f.write("piFlashNu                & " + str(self._zeroWeight['piFlashNu']) + " & " + str(self._NZWeight['piFlashNu']) + " &\\\\ \n")
        f.write("muonDecay                & " + str(self._zeroWeight['muonDecay']) + " & " + str(self._NZWeight['muonDecay']) + " &\\\\ \n")
        f.write("eProduction              & " + str(self._zeroWeight['eProduction']) + " & " + str(self._NZWeight['eProduction']) + " &\\\\ \n")
        f.write("$\\nu_{\\mu}$ production & " + str(self._zeroWeight['numuProduction']) + " & " + str(self._NZWeight['numuProduction']) + " &\\\\ \n")
        f.write("$\\nu_e$ production      & " + str(self._zeroWeight['nueProduction']) + " & " + str(self._NZWeight['nueProduction']) + " &\\\\ \n")
        f.write("$\\nu_{\\mu}$ detector   & " + str(self._zeroWeight['numuDetector']) + " & " + str(self._NZWeight['numuDetector']) + " &\\\\ \n")
        f.write("$\\nu_e$ detector        & " + str(self._zeroWeight['numuDetector']) + " & " + str(self._NZWeight['numuDetector']) + " &\\\\ \n")
        f.write("\\hline\n")
        f.write("\\end{tabular}  \\nl\n")

        f.close()
