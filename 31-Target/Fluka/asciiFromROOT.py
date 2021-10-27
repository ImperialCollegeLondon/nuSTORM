#!/usr/bin/env python3

"""
Python script to store the pion and muon tracking information in
ASCII format from ROOT files containing the (converted) Fluka output
from getPiMuDist.py. Each text line will contain the following variables:
"x y z px py pz PDGId TotE"
Weights are not written, but they should all be equal to unity anyway.
Units are: position = cm, momentum = GeV/c, energy = GeV.
BDSIM ASCII input beam distribution file info:
http://www.pp.rhul.ac.uk/bdsim/manual/model_control.html#userfile

"""

import os
import ROOT
import sys

def run(rootFileName, ntupleName = 'Data'):

    print('rootFileName is {0}'.format(rootFileName))

    # Equivalent text file
    textFileName = rootFileName.replace('.root', '.txt')
    print('Text file will be {0}'.format(textFileName))
    textFile = open(textFileName, 'w')

    # ROOT file and data
    rootFile = ROOT.TFile(rootFileName, 'read')
    data = rootFile.Get(ntupleName)

    # Loop over tree entries
    N = data.GetEntries()
    for i in range(N):

        if i%500000 == 0:
            print('Event {0}'.format(N - i))

        data.GetEntry(i)

        # Write out position, momentum, PDGId and total energy
        x = getattr(data, 'x')
        y = getattr(data, 'y')
        z = getattr(data, 'z')
        px = getattr(data, 'px')
        py = getattr(data, 'py')
        pz = getattr(data, 'pz')
        PDGId = getattr(data, 'PDGId')
        totE = getattr(data, 'totE')

        line = '{0:<.6e} {1:<.6e} {2:<.6e} {3:<.6e} {4:<.6e} {5:<.6e} ' \
                '{6:<} {7:<.6e}\n'.format(x, y, z, px, py, pz, PDGId, totE)
        textFile.write(line)

    textFile.close()
    rootFile.Close()
    # gzip the text file
    print('gzipping {0}'.format(textFileName))
    os.system('gzip {0}'.format(textFileName))


if __name__ == '__main__':

    # Base directory location for the ROOT files
    baseDir = 'fluka4-1.1/nuSTORM'
    rootFileName = 'nuSTORMTarget_26GeV_HPosAll/nuSTORMTarget_26GeV_HPos.root'
    NArg = len(sys.argv)

    if NArg > 1:
        rootFileName = sys.argv[1]

    run(rootFileName)
