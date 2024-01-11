#       singleTargetBin.py
#   Version 1.0                 Paul Kyberd
#
#   Create just one bin - based on Marvin's code which does all bins, but with integer means
#

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from iminuit import Minuit
from iminuit.cost import LeastSquares
import ROOT
import argparse

pStart = 1
pStop  = 10

# x = np.array([])
# y = np.array([])
# z = np.array([])
# xp = np.array([])
# yp = np.array([])
# ptot = np.array([])

parser = argparse.ArgumentParser()
parser.add_argument('--inFile', help='Select input file. Default: Fluka100GeVHPos.root', default='Fluka100GeVHPos.root')
parser.add_argument('--ntuple', help='Select ntuple name. Default: Data', default='Data')
parser.add_argument('--meanE', help='Set the mean Energy of the pions to extract', default=4)

#parser.add_argument('--p0', help='Select central pion momentum for histogram generation. Default: 5 [GeV].',default='5')
#parser.add_argument('--outFile', help='Select output file. Default: Scratch/target.root', default='Scratch/target.root')
parser.add_argument('--outFile', help='Select output file. Default: target (.root) will be appended', default='target')        # my mod
args = parser.parse_args()

print("=============  Histogram generation from ROOT file started  =============")

nuSIMPATH = os.getenv('nuSIMPATH')

print("------- nuSIMPATH is ", nuSIMPATH)

rootFileName = args.inFile
#rootOutFileName = os.path.join(nuSIMPATH, args.outFile)                    # mod Marvin's defaults
#   Mean Energy is 
meanE = float(args.meanE)
#   make outputfile name
#   Energy GeV
geV = int(meanE)
print("Integer part of mean Energy: ",geV)
#   fractional part
fracE = int((meanE*10 - geV*10)*10) 
print("fractional part of mean Energy: ",fracE)
#   Output file
rootOutFileName = args.outFile+str(geV)+str(fracE)+'.root'                  # into the target directory - labelled by mean Energy
ntupleName = args.ntuple

#   Mean Energy is 
meanE = float(args.meanE)

print("Input file: ",rootFileName)
print("Output file: ",rootOutFileName)
print("ntuple name: ",ntupleName)
print("mean Energy is : ",meanE)


rootFile = ROOT.TFile(rootFileName, 'read')
data = rootFile.Get(ntupleName)

n = data.GetEntries()

print("Input file opened...")
print("Entries in input file: ",n)

rootOutFile = ROOT.TFile(rootOutFileName, 'RECREATE', 'ROOT file with Histograms' )

print("Output file opened...")

pLow = []
pHigh = []

histsX = []
histsY = []
histsE = []
histsP = []
histsPx = []
histsPy = []
histsPz = []
histsXp = []
histsYp = []
histsXPS = []
histsYPS = []
histsXMS = []
histsYMS = []

#for j in range(pStart,pStop+1):

#   Create histograms - limits given by mean energy of the pions
#meanE = 3.5

pLow.append(meanE*0.9)
pHigh.append(meanE*1.1)
histsP.append(ROOT.TH1D('histP'+str(geV)+str(fracE)+'GeV', 'Momentum Histogram', 50, meanE*0.9, meanE*1.1))
histsE.append(ROOT.TH1D('histE'+str(geV)+str(fracE)+'GeV', 'Energy Histogram', 75, meanE*0.85, meanE*1.15))
histsX.append(ROOT.TH1D('histX'+str(geV)+str(fracE)+'GeV', 'X Histogram', 50, -0.1, 0.1))
histsY.append(ROOT.TH1D('histY'+str(geV)+str(fracE)+'GeV', 'Y Histogram', 50, -0.1, 0.1))
histsPx.append(ROOT.TH1D('histPx'+str(geV)+str(fracE)+'GeV', 'Px Histogram', 60, -0.6, 0.6))
histsPy.append(ROOT.TH1D('histPy'+str(geV)+str(fracE)+'GeV', 'Py Histogram', 60, -0.6, 0.6))
histsPz.append(ROOT.TH1D('histPz'+str(geV)+str(fracE)+'GeV', 'Pz Histogram', 75, meanE*0.80, meanE*1.1))
histsXp.append(ROOT.TH1D('histXp'+str(geV)+str(fracE)+'GeV', 'Xp Histogram', 100, -0.2, 0.2))
histsYp.append(ROOT.TH1D('histYp'+str(geV)+str(fracE)+'GeV', 'Yp Histogram', 100, -0.2, 0.2))
histsXPS.append(ROOT.TH2D('histXPS'+str(geV)+str(fracE)+'GeV', 'X PS Histogram', 200, -0.12, 0.12, 200, -0.35, 0.35))
histsYPS.append(ROOT.TH2D('histYPS'+str(geV)+str(fracE)+'GeV', 'Y PS Histogram', 200, -0.12, 0.12, 200, -0.35, 0.35))
histsXMS.append(ROOT.TH2D('histXMS'+str(geV)+str(fracE)+'GeV', 'X MS Histogram', 200, -0.12, 0.12, 200, -1.2, 1.2))
histsYMS.append(ROOT.TH2D('histYMS'+str(geV)+str(fracE)+'GeV', 'Y MS Histogram', 200, -0.12, 0.12, 200, -1.2, 1.2))

#histsP[0] = ROOT.TH1D('histP1GeV', 'Momentum Histogram', 20, 0.9, 1.1)
print("Histograms initialised... Starting reading data and filling histograms... ")

#   Loop over the input file
for i in range(n):

    data.GetEntry(i)

    x = getattr(data, 'x')
    y = getattr(data, 'y')
    #p = getattr(data, 'p')
    px = getattr(data, 'px')
    py = getattr(data, 'py')
    pz = getattr(data, 'pz')
    p = np.sqrt(px**2+py**2+pz**2)
    PDGId = getattr(data, 'PDGId')
    totE = getattr(data, 'totE')

    if PDGId == 211:
        for k in range(len(pLow)):
            if ((p>=pLow[k]) and (p<=pHigh[k])):

                rootOutFile.cd()

                histsP[k].Fill(p)
                histsE[k].Fill(totE)
                histsX[k].Fill(x/100.)
                histsY[k].Fill(y/100.)
                histsPx[k].Fill(px)
                histsPy[k].Fill(py)
                histsPz[k].Fill(pz)
                histsXp[k].Fill(px/pz)
                histsYp[k].Fill(py/pz)
                histsXPS[k].Fill(x/100.,px/pz)
                histsYPS[k].Fill(y/100.,py/pz)
                histsXMS[k].Fill(x/100.,px)
                histsYMS[k].Fill(y/100.,py)

        rootFile.cd()

    if (i < 10):
        print ("event number is ", i)
    elif ((i <100) and (i%10 ==0)):
        print ("event number is ", i)
    elif ((i <1000) and (i%100 ==0)):
        print ("event number is ", i)
    elif ((i <10000) and (i%1000 ==0)):
        print ("event number is ", i)
    elif ((i <100000) and (i%10000 ==0)):
        print ("event number is ", i)
    elif ((i <1000000) and (i%100000 ==0)):
        print ("event number is ", i)
    elif (i%1000000 ==0):
        print ("event number is ", i)

rootFile.Close()
rootOutFile.cd()

print("Data read, histograms filled & input file closed...")
print("Writing out histograms...")

for l in range(len(histsP)):
    histsP[l].Write()
    histsE[l].Write()
    histsX[l].Write()
    histsY[l].Write()
    histsPx[l].Write()
    histsPy[l].Write()
    histsPz[l].Write()
    histsXp[l].Write()
    histsYp[l].Write()
    histsXPS[l].Write()
    histsYPS[l].Write()
    histsXMS[l].Write()
    histsYMS[l].Write()

rootOutFile.Close()

print("=============  Histogram generation from ROOT file finished  =============")

# #! PLOTTING
# n, bins, patches = plt.hist(ptot, bins=50, color='y', range=(4.5,5.5))
# plt.xlabel('Momentum (GeV)')
# plt.ylabel('Entries')
# plt.title('Momentum distribution')
# plt.savefig('Scratch/psDistributions/0_target/ptot.pdf')
# plt.close()
#
# n, bins, patches = plt.hist(x, bins=50, color='y', range=(-10.,10.))  #, range=(-0.16,0.16)   , range=(-0.3,0.3)
# plt.xlabel('x [m]')
# plt.ylabel('Entries')
# plt.title('x distribution')
# plt.savefig('Scratch/psDistributions/0_target/x.pdf')
# plt.close()
#
# n, bins, patches = plt.hist(y, bins=50, color='y', range=(-10.,10.))  #, range=(-0.16,0.16)   , range=(-0.3,0.3)
# plt.xlabel('y [m]')
# plt.ylabel('Entries')
# plt.title('y distribution')
# plt.savefig('Scratch/psDistributions/0_target/y.pdf')
# plt.close()
#
# n, bins, patches = plt.hist(z, bins=50, color='y', range=(-0.01,0.01))
# plt.xlabel('z [m]')
# plt.ylabel('Entries')
# plt.title('z distribution')
# plt.savefig('Scratch/psDistributions/0_target/z.pdf')
# plt.close()
#
# n, bins, patches = plt.hist(xp, bins=50, color='y', range=(-0.35,0.35)) #, range=(-0.0075,0.0075)   , range=(-0.01,0.01)
# plt.xlabel('x^prime')
# plt.ylabel('Entries')
# plt.title('x^prime distribution')
# plt.savefig('Scratch/psDistributions/0_target/xp.pdf')
# plt.close()
#
# n, bins, patches = plt.hist(yp, bins=50, color='y', range=(-0.35,0.35)) #, range=(-0.0075,0.0075)   , range=(-0.01,0.01)
# plt.xlabel('y^prime')
# plt.ylabel('Entries')
# plt.title('y^prime distribution')
# plt.savefig('Scratch/psDistributions/0_target/yp.pdf')
# plt.close()
#
# pltX = plt.hist2d(x,xp, bins=200,norm=colors.LogNorm(),range=[[-10.,10.],[-.35,.35]]) #,range=[[-0.5,0.5],[-.05,.05]]
# plt.colorbar()
# plt.xlabel('x [m]')
# plt.ylabel('x^prime')
# plt.title('phase space distribution in x')
# #plt.grid(alpha=0.)
# #plt.show()
# plt.savefig('Scratch/psDistributions/0_target/xPS.pdf')
# plt.close()
#
# pltY = plt.hist2d(y,yp, bins=200,norm=colors.LogNorm(),range=[[-10.,10.],[-.35,.35]]) #, range=[[-0.5,0.5],[-.05,.05]]
# plt.colorbar()
# plt.xlabel('y [m]')
# plt.ylabel('y^prime')
# plt.title('phase space distribution in y')
# #plt.grid(False)
# #plt.show()
# plt.savefig('Scratch/psDistributions/0_target/yPS.pdf')
# plt.close()