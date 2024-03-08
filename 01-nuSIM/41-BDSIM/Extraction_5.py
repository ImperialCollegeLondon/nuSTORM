import sys
import pybdsim
import os.path as _path
import numpy as _np

def _ExtractSamplerCoords(data, samplername):
    """ Extract sampler coordinates."""
    if samplername != "Primary":
        # add . to the sampler name to match branch names from file
        if samplername[-1] != ".":
            samplername += "."
        # check branch exists
        allSamplers = data.GetSamplerNames()
        if not samplername in allSamplers:
            print("Sampler " + samplername + " not found in inputfile. Terminating...")
            sys.exit(1)

    sampler = pybdsim.Data.SamplerData(data, samplername)

    # get particle coords and filter out secondaries
    x   = sampler.data['x']
    y   = sampler.data['y']
    xp = sampler.data['xp']
    yp  = sampler.data['yp']
    tof = sampler.data['T']
    E   = sampler.data['energy']
    pid = sampler.data['partID']

    return x,xp,y,yp,tof,E,pid


def BdsimSampler2BdsimUserfile(inputfile, outfile, samplername):
    """
    Takes .root file generated from a BDSIM run and creates
    a BDSIM userFile file from the sampler particle tree.
    inputfile   - <str> root format output from BDSIM run
    outfile     - <str> filename for the inrays file
    samplername - <str> sampler name in BDSIM root file
    Writes the following columns to file:
      x[m] xp[rad] y[m] yp[rad] t[ns] E[GeV] PID[]
    E is the total particle energy.
    The t column is the time in the given sampler minus the mean time for that sampler.
    If not mean subtracted, the particles may be significantly offset from the primary position.
    """
    if not (outfile[-4:] == ".dat"):
        outfile = outfile + ".dat"

    if isinstance(inputfile, str):
        if not _path.isfile(inputfile):
            raise IOError("file \"{}\" not found!".format(inputfile))
        else:
            print("Loading input file: ", inputfile)
            data = pybdsim.Data.Load(inputfile)

    x,xp,y,yp,tof,E,pid = _ExtractSamplerCoords(data, samplername)

    # subtract mean time as non-primary sampler will be at a finite T in the lattice - should be centred around 0.
    if samplername != "Primary":
        meanT = _np.mean(tof)
        tof = tof - meanT
    nentries = len(x)

    outfile = open(outfile, 'w')
    for n in range(0, nentries):  # n denotes a given particle
        s =  ' ' + repr(x[n])
        s += ' ' + repr(y[n])
        s += ' ' + repr(xp[n])
        s += ' ' + repr(yp[n])
        s += ' ' + repr(tof[n])
        s += ' ' + repr(E[n])
        s += ' ' + repr(pid[n])
        s += '\n'
        outfile.writelines(s)

    outfile.close()

BdsimSampler2BdsimUserfile('/vols/ccap/users/tta20/tta20/testRohan.root', 'testRohan.dat', samplername='d35_5')
# BdsimSampler2BdsimUserfile('Complete_Line_5GeV_nQuads.root','particlesCoords_5GeV_nQuads.dat', samplername='d35_5')

