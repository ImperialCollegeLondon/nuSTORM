#!/usr/bin/env python
#       histoManager.py                                    Version 2.0
#
#     Paul Kyberd                                               29 October 2021
#
#  Version 2.0                                                  29 October 2021
#  Update to python3
#
#  Version 1.0                                                  25 March 2013
# Attempt to get the booking and filling of histograms rather simpler to manage
# ie go back to hbook
#
import ROOT
from testUnit import testUnit
from pathlib import Path

class histo():
    """
    Histogram class - facade
    """

    def __init__(self):
        self.Version = 1.0
        self.hists = 0

    def status(self):
        return False

    def version(self):
        return self.Version


    def book(self, title, bins, lower, upper):
#
#    babar page shows the call as
#  TH1D *h1 = new TH1d("h1", "title", bins, lower, upper           - do the *h1 and h1 have to be the same
#  TH2D ... bins, lower,upper, bins,lower , upper
#  TH3D ... bins, lower,upper, bins,lower,upper, bins,lower,upper
#
#         Start by just doing a simple 1D

        self.Title = title
        name = title                                  # name must be unique .. put it equal to title
        self.histVar = ROOT.TH1D(name, title, bins, lower, upper)

    def Fill(self, value):
        self.histVar.Fill(value)

    def output(self):
        canvas = ROOT.TCanvas(self.Title, self.Title)
        self.histVar.Draw()
        canvas.Draw()
        canvas.Print(self.Title+".png")

class histoManager:
    """
    Create and manage histograms
    """


    def __init__(self):
        self.Version = 1.0
        self.histNames = []
        self.hists=[]
        self.histParams=[]
        self.action = ""

    def toString(self):

        retString = "Histogram manager " + str(self.Version) + "\nDebug Status         " + str(self.status())
        retString = retString + "\n" +  str(self.hists) + "\n" + str(self.histParams)

        return retString

    def book(self, title, bins, lower, upper):


#  Create a histo
#
#    babar page shows the call as
#  TH1D *h1 = new TH1d("h1", "title", bins, lower, upper           - do the *h1 and h1 have to be the same
#  TH2D ... bins, lower,upper, bins,lower , upper
#  TH3D ... bins, lower,upper, bins,lower,upper, bins,lower,upper
#
#         Start by just doing a simple 1D

        self.Title = title
        name = title                                  # name must be unique .. put it equal to title
        self.histVar = ROOT.TH1D(name, title, bins, lower, upper)
        self.hists.append(self.histVar)
        self.histParams.append("")
        return self.histVar

    def setParams(self,params):

        self.histParams[-1] = params

    def book2(self, title, bins1, lower1, upper1, bins2, lower2, upper2):

        self.Title = title
        name = title                                  # name must be unique .. put it equal to title
        self.histVar = ROOT.TH2D(name, title, bins1, lower1, upper1, bins2, lower2, upper2 )
        self.hists.append(self.histVar)
        self.histParams.append("")
        return self.histVar

    def histdo(self,path):
# Output all histos to screen and file
        hPnt = 0
        for i in range(len(self.hists)):
            hPnt = hPnt + 1
            hCurr = self.hists[i]
            pCurr = self.histParams[i]
            title = hCurr.GetTitle()
            if pCurr == "m":
                self.action = "m"
                hPnt = 1
                canvas = ROOT.TCanvas(title, title)
                canvas.Divide(2,2)
                canvas.cd(hPnt)
                canvas.Update()
                hCurr.Draw()

            if self.action == "m":
                canvas.cd(hPnt)
                canvas.Update()
                hCurr.Draw()
                if hPnt == 4:
                    canvas.Draw()
                    canvas.Print(path + "/plots/" + title + ".pdf")
                    hPnt = 0
                    self.action = ""
            else:
                canvas = ROOT.TCanvas(title, title)
                hCurr.Draw()
                canvas.Draw()
                canvas.Print(path + "/plots/" + title + ".pdf")

# create a tex file with all the plots

    def texCreate(self, fileName):
        texHline = "\\hline\n"
        slash = "\\"
#  Outout the information as a document
#        plotDir = Path.cwd()
#        plotFile = str(Path.cwd())+"/plots.tex"
        plotFile = fileName
        print (plotFile)
        f = open(plotFile, "w")
        f.write(slash + "documentclass{article}\n")
        f.write(slash + "usepackage{graphicx}\n")
        f.write("\n")
        f.write(slash + "begin{document}\n")
        f.write(slash + "title{Normalisation plots}\n")
        f.write("\n")

        hPnt = 0
        curHead = ""
        for i in range(len(self.hists)):
            hPnt = hPnt + 1
            hCurr = self.hists[i]
            pCurr = self.histParams[i]
            title = hCurr.GetTitle()
            end = title.find(":")
            newHead = title[0:end]
            print ("newHead is ", newHead)
            print ("title is ", title)

            if (curHead != newHead):
                curHead = newHead
                f.write("\n" + slash + "hfill" + slash + "newline\n")
                f.write("\n" + slash + "subsection*{" + curHead + "}\n" + slash + "noindent\n")
                hPnt = 1

            f.write(slash + "includegraphics[scale=0.40]{plots/" + title +  ".pdf}\n")
            if (hPnt%2 == 0): f.write("\n")

        f.write("\n")



        f.write(slash + "end{document}")
        f.close()



    def histOutRoot(self, fileName):
        self.outfile = ROOT.TFile( fileName, 'RECREATE', 'ROOT file with Histograms' )
        for i in range(len(self.hists)):
            hCurr = self.hists[i]
            hCurr.Write()


    def status(self):
        return True

    def version(self):
        return self.Version




if __name__ == "__main__":

# call the unit testing creation
    t = testUnit()

    hm = histoManager()
    t.announce(hm)
    t.fTry(["Version test",hm.version(),1.0])
    t.fTry(["Status test",hm.status(),True])
    stringComp = "Histogram manager 1.0\nDebug Status         True\n[]\n[]"
    t.fTry(["toString test - (empty)",hm.toString(),stringComp])

#    print (hm.toString())
#    print (stringComp)

    hTitle = "Test histogram 1"
    hBins  = 100
    hLower = -1.0
    hUpper = 1.0
    h1 = hm.book(hTitle, hBins, hLower, hUpper)
    hm.setParams("m")

    h3 = hm.book("Test histogram 3", 100, -3.0, 3.0)
    h4 = hm.book("Test histogram 4", 80, -4.0, 4.0)
    h5 = hm.book("Test histogram 5", 80, -4.0, 4.0)
    h6 = hm.book("Test histogram 6", 80, -4.0, 4.0)

    h3.Fill(0.3)
    h3.Fill(0.3)
    h4.Fill(0.4)
    h4.Fill(-0.5)
    h5.Fill(-0.5)
    h5.Fill(-0.5)
    h5.Fill(-0.5)
    h6.Fill(0.6)
    h6.Fill(0.6)
    h6.Fill(0.6)
    h6.Fill(0.6)

    print (hm.toString())
    hm.histdo()

    t.summary()
