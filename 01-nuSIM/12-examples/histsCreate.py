"""
Model for calculating normalised numbers
========================================

@date: 05 Jan 2022
@author: Paul Kyberd
@version: 2.0

@description: read in a the dictionary with the histogram values in from an external file

"""

from pathlib import Path
import math
import json
import numpy as np

class histsCreate():


    def __init__(self, hM, plotsDict):
        self._locsStart = []
        self._locs=[]
        self._hists=[]
        self._hm = hM
        self._nHists = 10
        self._count = 0
        self._targetWeight = 0
        self._targetNoWeight = 0
        self._t0 = 0
        self._NZWeight = {"target":0, "productionStraight": 0, "prodStraightEnd": 0, "pionDecay": 0, "muonProduction": 0,"piFlashNu": 0,
                "muonDecay": 0, "eProduction": 0, "numuProduction":0, "nueProduction": 0, "numuDetector": 0, "nueDetector": 0}
        self._zeroWeight = {"target":0, "productionStraight": 0, "prodStraightEnd": 0, "pionDecay": 0, "muonProduction": 0,"piFlashNu": 0,
                "muonDecay": 0, "eProduction": 0, "numuProduction":0, "nueProduction": 0, "numuDetector": 0, "nueDetector": 0}

#  read in the dictionaries
        with open(plotsDict) as pD:
            self._plotsInfo = json.load(pD)

    def __repr__(self):
        return "create a set of histograms for a particle at a location in the event History"

    def histAdd(self, eventType):
        self._locsStart.append(len(self._hists))
        self._locs.append(eventType)
        print ("locStart is ", self._locsStart)
        print ("locs ", self._locs)
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

#   r v rx
        hTitle = eventType + ":x v px"
        hLow1 = self._plotsInfo["xLower"][eventType]
        hUp1 = self._plotsInfo["xHigher"][eventType]
        hLow2 = self._plotsInfo["pxLower"][eventType]
        hUp2 = self._plotsInfo["pxHigher"][eventType]
        self._hists.append(self._hm.book2(hTitle, 100, hLow1, hUp1, 100, hLow2, hUp2))
        hTitle = eventType + ":y v py"
        hLow1 = self._plotsInfo["yLower"][eventType]
        hUp1 = self._plotsInfo["yHigher"][eventType]
        hLow2 = self._plotsInfo["pyLower"][eventType]
        hUp2 = self._plotsInfo["pyHigher"][eventType]
        self._hists.append(self._hm.book2(hTitle, 100, hLow1, hUp1, 100, hLow2, hUp2))

#   phase space

        hTitle = eventType + ":x v xp"
        hLow1 = self._plotsInfo["xLower"][eventType]
        hUp1 = self._plotsInfo["xHigher"][eventType]
#       hLow2 = self._plotsInfo["xpLower"][eventType]
#       hUp2 = self._plotsInfo["xpHigher"][eventType]
        hLow2 = -0.05
        hUp2 = 0.05
        print ("xpLower is ", hLow2, "    xpHigher is ", hUp2)
        self._hists.append(self._hm.book2(hTitle, 100, hLow1, hUp1, 100, hLow2, hUp2))
        hTitle = eventType + ":y v yp"
        hLow1 = self._plotsInfo["yLower"][eventType]
        hUp1 = self._plotsInfo["yHigher"][eventType]
        hLow2 = self._plotsInfo["ypLower"][eventType]
        hUp2 = self._plotsInfo["ypHigher"][eventType]
        self._hists.append(self._hm.book2(hTitle, 100, hLow1, hUp1, 100, hLow2, hUp2))

#   correlations
        hTitle = eventType + ":x v y"
        hLow1 = self._plotsInfo["xLower"][eventType]
        hUp1 = self._plotsInfo["xHigher"][eventType]
        hLow2 = self._plotsInfo["yLower"][eventType]
        hUp2 = self._plotsInfo["yHigher"][eventType]
        self._hists.append(self._hm.book2(hTitle, 100, hLow1, hUp1, 100, hLow2, hUp2))

        hTitle = eventType + ":xp v yp"
        hLow1 = self._plotsInfo["xpLower"][eventType]          # don't check here, should have alread been caught
        hUp1 = self._plotsInfo["xpHigher"][eventType]
        hLow2 = self._plotsInfo["ypLower"][eventType]
        hUp2 = self._plotsInfo["ypHigher"][eventType]
        self._hists.append(self._hm.book2(hTitle, 100, hLow1, hUp1, 100, hLow2, hUp2))


#   meant to show the bunch structure
        hLower = self._plotsInfo["tBunchLower"][eventType]
        hUpper = self._plotsInfo["tBunchHigher"][eventType]
        hTitle = eventType + ":t from bunch start"
        if eventType == "muonDecay" or eventType == "eProduction":
            hBins = 50
        elif eventType == "target":
            hbins = 70
            hUpper = 70.
        else:
            hBins = 100
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))

#   meant to show the spill structure
        hLower = self._plotsInfo["tSpillLower"][eventType]
        hUpper = self._plotsInfo["tSpillHigher"][eventType]
        hTitle = eventType + ":t for spill"
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))
        hLower = self._plotsInfo["TOFLower"][eventType]
        hUpper = self._plotsInfo["TOFHigher"][eventType]
        hTitle = eventType + ":Time of Flight"
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
        print ("hLower for pz is ", hLower, "   hUpper for pz is ", hUpper)
        print ("hTtile is ", hTitle)
        hist = self._hm.book(hTitle, hBins, hLower, hUpper)
        print ("hist is ", hist)
        self._hists.append(hist)
#  Energy plot
        hTitle = eventType + ":p"
        try:
            hLower = self._plotsInfo["pLower"][eventType]
            hUpper = self._plotsInfo["pHigher"][eventType]
        except:
            print ("error is ", OSError)
            exit("Change plotsXXX.dict ... probably need to replace eLower and eHigher by pLower and pHigher")
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))

        hTitle = eventType + ":p_v_t"
        hLow1 = self._plotsInfo["pLower"][eventType]
        hUp1 = self._plotsInfo["pHigher"][eventType]
        hLow2 = self._plotsInfo["TOFLower"][eventType]
        hUp2 = self._plotsInfo["TOFHigher"][eventType]
        self._hists.append(self._hm.book2(hTitle, 50, hLow1, hUp1, 50, hLow2, hUp2))

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
        p = np.sqrt(px*px + py*py +pz*pz)
        wt = particle.weight()
        s = particle.s()
        t = particle.t()

#  if the location is target get the time as the start time
        if (location == 'target'):
            self._t0 = t

        if (wt !=0):
            n = self._NZWeight.get(location)
            self._NZWeight[location] = n + 1
        else:
            n = self._zeroWeight.get(location)
            self._zeroWeight[location] = n + 1
#   book order ... fill must be same
#0;#1; 2;      3; 4;      5;      6;      7;      8;     9;      10;     11;     12;  13; 14; 15; 16;17;    18;   
#x; y; z; weight; s; x v px; y v py; x v xp; y v yp; x v y; xp v yp; tbunch; tSpill; TOF; px; py; pz; p; p v t;

        hoffset = 0
        if (wt > 0.0):
            self._hists[hPnt+hoffset].Fill(x)       # 0
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(y)       # 1
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(z)       # 2
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(wt)      # 3
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(s)       # 4
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(x,px)    # 5
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(y,py)    # 6
            hoffset = hoffset + 1
            xp = px/pz
            yp = py/pz
            self._hists[hPnt+hoffset].Fill(x,xp)    # 7
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(y,yp)    # 8
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(x,y)     # 9
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(xp,yp)   # 10
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(t)       # 11
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(t-self._t0)  # 12
            hoffset = hoffset + 1
            # TOF                                   # 13            
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(px)      # 14
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(py)      # 15
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(pz)      # 16
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(p)      # 17
            hoffset = hoffset + 1
            self._hists[hPnt+hoffset].Fill(pz,t)    # 18

            if (location == 'productionStraight'):
# plots with non-zero weight - for production straight
                 self._count = self._count + 1

            Enu = math.sqrt(px*px + py*py + pz*pz)
            if ((location == 'numuDetector') or (location == "nueDetector")):
                if ((abs(x) < 2.5) and (abs(y) < 2.5)):
                        self._hists[hPnt+13].Fill(pz)
                        self._hists[hPnt+14].Fill(Enu)
            else:
                self._hists[hPnt+13].Fill(pz)
                self._hists[hPnt+14].Fill(Enu)



    def summary(self, fileName):

        texHline = "\\hline\n"
        slash = "\\"
        print(" Summary of results")
        print (self._zeroWeight)
        print (self._NZWeight)
#  Outout the information as a latex table
#        sumDir = Path.cwd()
#        sumFile = str(Path.cwd())+"/summary.tex"
        sumFile = fileName
        print (sumFile)
        f = open(sumFile, "w")
        f.write("\\begin{tabular}{|  l | l | l  | p{9.0cm} | }\n" )
        f.write(texHline )
        f.write("Location   &  zero weight  &  Full weight & comments \\\\ \n"  )
        f.write(texHline )
        f.write("target                   & " + str(self._zeroWeight['target']) + " & " + str(self._NZWeight['target']) + " &\\\\ \n")
        f.write("productionStraight       & " + str(self._zeroWeight['productionStraight']) + " & " + str(self._NZWeight['productionStraight']) + " &\\\\ \n")
        f.write("prodStraightEnd          & " + str(self._zeroWeight['prodStraightEnd']) + " & " + str(self._NZWeight['prodStraightEnd']) + " &\\\\ \n")
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
