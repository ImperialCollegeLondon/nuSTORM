#       singleTargetBin2.py
#   Version 1.0                 Paul Kyberd
#
#   Modified version of singleTargetBin ... to read in the .txt version of Paul Jurj's file
#   The values are:
#   The x position "x" (cm)
#   The y position "y" (cm)
#   The z position "z" (cm)
#   The x momentum "px" (GeV/c)
#   The y momentum "py" (GeV/c)
#   The z momentum "pz" (GeV/c)
#   The x angle (rad)
#   The y angle (rad)
#   Total energy "totE" (GeV)
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
parser.add_argument('--inFile', help='Select input file. Default: selected_pions_5GeV_1.0mm.txt', default='selected_pions_5GeV_1.0mm.txt')
parser.add_argument('--ntuple', help='Select ntuple name. Default: Data', default='Data')
parser.add_argument('--meanE', help='Set the mean Energy of the pions to extract', default=5)

#parser.add_argument('--p0', help='Select central pion momentum for histogram generation. Default: 5 [GeV].',default='5')
#parser.add_argument('--outFile', help='Select output file. Default: Scratch/target.root', default='Scratch/target.root')
parser.add_argument('--outFile', help='Select output file. Default: target (.root) will be appended', default='target')        # my mod
args = parser.parse_args()

print("=============  Histogram generation from .txt file started (Paul Jurj format) =============")

nuSIMPATH = os.getenv('nuSIMPATH')

print("------- nuSIMPATH is ", nuSIMPATH)

txtFileName = args.inFile
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

print("Input file: ",txtFileName)
print("Output file: ",rootOutFileName)
print("ntuple name: ",ntupleName)
print("mean Energy is : ",meanE)

rootOutFile = ROOT.TFile(rootOutFileName, 'RECREATE', 'ROOT file with Histograms' )

print("Output file opened...")

pLow = []
pHigh = []

histsX = []
histsY = []
histsE = []
histsPx = []
histsPy = []
histsPz = []
histsXp = []
histsYp = []
histsXPS = []
histsYPS = []
histsXMS = []
histsYMS = []

pLow.append(meanE*0.9)
pHigh.append(meanE*1.1)
histsP = ROOT.TH1D('histP'+str(geV)+str(fracE)+'GeV', 'Momentum Histogram', 50, meanE*0.9, meanE*1.1)
histsE = ROOT.TH1D('histE'+str(geV)+str(fracE)+'GeV', 'Energy Histogram', 75, meanE*0.85, meanE*1.15)
histsX = ROOT.TH1D('histX'+str(geV)+str(fracE)+'GeV', 'X Histogram', 50, -0.1, 0.1)
histsY = ROOT.TH1D('histY'+str(geV)+str(fracE)+'GeV', 'Y Histogram', 50, -0.1, 0.1)
histsPx = ROOT.TH1D('histPx'+str(geV)+str(fracE)+'GeV', 'Px Histogram', 60, -0.6, 0.6)
histsPy = ROOT.TH1D('histPy'+str(geV)+str(fracE)+'GeV', 'Py Histogram', 60, -0.6, 0.6)
histsPz = ROOT.TH1D('histPz'+str(geV)+str(fracE)+'GeV', 'Pz Histogram', 75, meanE*0.80, meanE*1.2)
histsXp = ROOT.TH1D('histXp'+str(geV)+str(fracE)+'GeV', 'Xp Histogram', 100, -0.2, 0.2)
histsYp = ROOT.TH1D('histYp'+str(geV)+str(fracE)+'GeV', 'Yp Histogram', 100, -0.2, 0.2)
histsXPS = ROOT.TH2D('histXPS'+str(geV)+str(fracE)+'GeV', 'X PS Histogram', 200, -0.12, 0.12, 200, -0.35, 0.35)
histsYPS = ROOT.TH2D('histYPS'+str(geV)+str(fracE)+'GeV', 'Y PS Histogram', 200, -0.12, 0.12, 200, -0.35, 0.35)
histsXMS = ROOT.TH2D('histXMS'+str(geV)+str(fracE)+'GeV', 'X MS Histogram', 200, -0.12, 0.12, 200, -1.2, 1.2)
histsYMS = ROOT.TH2D('histYMS'+str(geV)+str(fracE)+'GeV', 'Y MS Histogram', 200, -0.12, 0.12, 200, -1.2, 1.2)

print("Histograms initialised... Starting reading data and filling histograms... ")


f = open(txtFileName, "r")
print("Input file opened...")

#   Loop over the input file

i = 0
for line in f.readlines():


    data = line.split()

    x = float(data[0])
    y = float(data[1])
    z = float(data[2])
    px = float(data[3])
    py = float(data[4])
    pz = float(data[5])
    xpr = float(data[6])
    ypr = float(data[7])
    totE = float(data[8])
    p = np.sqrt(px**2+py**2+pz**2)

    PDGId = 211
    rootOutFile.cd()


    histsP.Fill(p)
    histsE.Fill(totE)
    histsX.Fill(x)
    histsY.Fill(y)
    histsPx.Fill(px)
    histsPy.Fill(py)
    histsPz.Fill(pz)
    histsXp.Fill(xpr)
    histsYp.Fill(ypr)
    histsXPS.Fill(x,xpr)
    histsYPS.Fill(y,ypr)
    histsXMS.Fill(x,px)
    histsYMS.Fill(y,py)

    i = i + 1
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


f.close()


rootOutFile.cd()

print("Data read, histograms filled & input file closed...")
print("Writing out histograms...")


histsP.Write()
histsE.Write()
histsX.Write()
histsY.Write()
histsPx.Write()
histsPy.Write()
histsPz.Write()
histsXp.Write()
histsYp.Write()
histsXPS.Write()
histsYPS.Write()
histsXMS.Write()
histsYMS.Write()

rootOutFile.Close()

print("=============  Histogram generation from text file finished  =============")

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