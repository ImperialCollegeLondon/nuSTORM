from pathlib import Path

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

        self._xLower = {"target":-0.5, "productionStraight":-1.0, "pionDecay":-2.0, "muonProduction":-1.0,"piFlashNu":-2.0,
                "muonDecay":-1.0, "eProduction":-1.0, "numuProduction":-1.0, "nueProduction":-1.0, "numuDetector":-50.0, "nueDetector":-10.0}
        self._xHigher = {"target": 0.5, "productionStraight": 1.0, "pionDecay": 100.0, "muonProduction": 1.0,"piFlashNu": 100.0,
                "muonDecay": 1.0, "eProduction": 1.0, "numuProduction":0, "nueProduction": 1.0, "numuDetector": 150.0, "nueDetector": 10.0}
        self._yLower = {"target":-0.5, "productionStraight":-1.0, "pionDecay":-1.0, "muonProduction":-1.0,"piFlashNu":-1.0,
                "muonDecay": -1.0, "eProduction":-1.0, "numuProduction":0, "nueProduction":-1.0, "numuDetector":-100.0, "nueDetector":-10.0}
        self._yHigher = {"target":0.5, "productionStraight": 1.0, "pionDecay": 1.0, "muonProduction": 1.0,"piFlashNu": 1.0,
                "muonDecay": 1.0, "eProduction": 1.0, "numuProduction":0, "nueProduction": 1.0, "numuDetector": 100.0, "nueDetector": 10.0}
        self._zLower = {"target": -60.0, "productionStraight": 0.0, "pionDecay": -60.0, "muonProduction": -100.0,"piFlashNu": -80.0,
                "muonDecay": -100.0, "eProduction": -100.0, "numuProduction": -100.0, "nueProduction": -100.0, "numuDetector": 0.0, "nueDetector":-1.0}
        self._zHigher = {"target":-40.0, "productionStraight": 100.0, "pionDecay": 200.0, "muonProduction": 500.0,"piFlashNu": 50.0,
                "muonDecay": 400.0, "eProduction": 400.0, "numuProduction": 400.0, "nueProduction": 400.0, "numuDetector": 500.0, "nueDetector": 500.0}
        self._pxLower = {"target":-0.2, "productionStraight":-0.2, "pionDecay":-1.0, "muonProduction":-1.0,"piFlashNu":-10.0,
                "muonDecay":-1.0, "eProduction":-1.0, "numuProduction":-1.0, "nueProduction":-1.0, "numuDetector":-10.0, "nueDetector":-10.0}
        self._pxHigher = {"target": 0.2, "productionStraight": 0.2, "pionDecay": 1.0, "muonProduction": 1.0,"piFlashNu": -20.0,
                "muonDecay": 1.0, "eProduction": 1.0, "numuProduction":0, "nueProduction": 1.0, "numuDetector": 10.0, "nueDetector": 10.0}
        self._pyLower = {"target":-0.2, "productionStraight":-0.2, "pionDecay":-1.0, "muonProduction":-1.0,"piFlashNu": 20.0,
                "muonDecay":-1.0, "eProduction":-1.0, "numuProduction":-1.0, "nueProduction":-1.0, "numuDetector":-10.0, "nueDetector":-10.0}
        self._pyHigher = {"target": 0.2, "productionStraight": 0.2, "pionDecay": 1.0, "muonProduction": 1.0,"piFlashNu": 1.0,
                "muonDecay": 1.0, "eProduction": 1.0, "numuProduction":0, "nueProduction": 1.0, "numuDetector": 10.0, "nueDetector": 10.0}
        self._pzLower = {"target": 0.0, "productionStraight": 0.0, "pionDecay":-1.0, "muonProduction":-1.0,"piFlashNu":-1.0,
                "muonDecay":-1.0, "eProduction":-6.0, "numuProduction":-6.0, "nueProduction":-6.0, "numuDetector":-10.0, "nueDetector":-10.0}
        self._pzHigher = {"target": 7.0, "productionStraight": 7.0, "pionDecay": 7.0, "muonProduction": 7.0,"piFlashNu": 7.0,
                "muonDecay": 7.0, "eProduction": 6.0, "numuProduction": 6.0, "nueProduction": 6.0, "numuDetector": 10.0, "nueDetector": 10.0}
        self._tLower = {"target": 0.0, "productionStraight": 0.0, "pionDecay": 0.0, "muonProduction": 0.0,"piFlashNu": 0.0,
                "muonDecay": 0.0, "eProduction": 0.0, "numuProduction":0, "nueProduction": 0.0, "numuDetector": 0.0, "nueDetector": 0.0}
        self._tHigher = {"target": 500.0, "productionStraight": 500.0, "pionDecay": 1000.0, "muonProduction": 500.0,"piFlashNu": 60.0,
                "muonDecay": 60000.0, "eProduction": 60000.0, "numuProduction": 60000.0, "nueProduction": 60000.0, "numuDetector": 60000.0, 
                "nueDetector": 60000.0}
        self._sLower = {"target":-1.0, "productionStraight": -1.0, "pionDecay": -1.0, "muonProduction": -1.0,"piFlashNu": -10.0,
                "muonDecay": 0.0, "eProduction": 0.0, "numuProduction":0.0, "nueProduction": 0.0, "numuDetector": 0.0, "nueDetector": 0.0}
        self._sHigher = {"target":400.0, "productionStraight": 400.0, "pionDecay": 1000.0, "muonProduction": 400.0,"piFlashNu": 70.0,
                "muonDecay": 60000.0, "eProduction": 60000.0, "numuProduction":40000.0, "nueProduction": 40000.0, 
                "numuDetector": 80000.0, "nueDetector": 80000.0}
        self._bins = {"target": 100, "productionStraight": 100, "pionDecay": 100, "muonProduction": 100,"piFlashNu": 100,
                "muonDecay": 50, "eProduction": 1000, "numuProduction":40000.0, "nueProduction": 40000.0, 
                "numuDetector": 40000.0, "nueDetector": 40000.0}


    def __repr__(self):
        return "create a set of histograms for a particle at a location in the event History"

    def histAdd(self, eventType):
        self._locsStart.append(len(self._hists))
        self._locs.append(eventType)
#        print ("locStart is ", self._locsStart)
#        print ("locs ", self._locs)
        hTitle = eventType + ":x" 
        hBins  = 100
        hLower = self._xLower[eventType]
        hUpper = self._xHigher[eventType]
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))

        hLower = self._yLower[eventType]
        hUpper = self._yHigher[eventType]
        hTitle = eventType + ":y" 
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))
        hLower = self._zLower[eventType]
        hUpper = self._zHigher[eventType]
        hTitle = eventType + ":z" 
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))
        hTitle = eventType + ":weight" 
        hLower = -10.0
        hUpper = 100.0
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))
        hLower = self._sLower[eventType]
        hUpper = self._sHigher[eventType]
        hTitle = eventType + ":s"
        if eventType == "muonDecay" or eventType == "eProduction":
            hBins = 50
        else:
            hBins = 100
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))
        hLower = self._tLower[eventType]
        hUpper = self._tHigher[eventType]
        hTitle = eventType + ":t" 
        if eventType == "muonDecay" or eventType == "eProduction":
            hBins = 50
        else:
            hBins = 100
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))
        hBins = 100
        hTitle = eventType + ":px" 
        hLower = self._pxLower[eventType]
        hUpper = self._pxHigher[eventType]
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))
        hTitle = eventType + ":py" 
        hLower = self._pyLower[eventType]
        hUpper = self._pyHigher[eventType]
        self._hists.append(self._hm.book(hTitle, hBins, hLower, hUpper))
        hTitle = eventType + ":pz" 
        hLower = self._pzLower[eventType]
        hUpper = self._pzHigher[eventType]
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
            self._hists[hPnt+8].Fill(pz)
            self._hists[hPnt+9].Fill(pz,t)

            if (location == 'productionStraight'):
# plots with non-zero weight - for production straight
                 self._count = self._count + 1

            if (location == "pionDecay"):
                if (s < 230.0): 
                    self._hists[hPnt+0].Fill(x)
                    self._hists[hPnt+1].Fill(y)
                    self._hists[hPnt+2].Fill(z)
            else:
                self._hists[hPnt+0].Fill(x)
                self._hists[hPnt+1].Fill(y)
                self._hists[hPnt+2].Fill(z)

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
