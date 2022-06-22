#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class RandomGenerator
=====================

  Class to generate random numbers according to a histogram.

  Class attributes:
  -----------------
  __instance : Set on creation of first (and only) instance.

  Instance attributes:
  --------------------
  _rootfilename = Filename, including path, to ROOT file containing
                  histgrams that should be used to generate random
                  numbers.
  _hist         = 1D histogram containing e.g. absolute momentum distribution
  _hist2Dx      = 2D histogram containing phase space distribution in x
  _hist2Dy      = 2D histogram containing phase space distribution in y

  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __new__ : Creates singleton class and prints version, PDG
                reference, and values of constants used.
      __repr__: One liner with call.
      __str__ : Dump of contents

  Get/set methods:
      CdVrsn()       : Returns code version number.
      getRootFilename: Returns filename of ROOT file, including path.

  Simulation methods:
      getRandom    : Generates one random number  according to a 1D histogram
      getRandom2Dx : Generates two random numbers according to a 2D histogram (in x)
      getRandom2Dy : Generates two random numbers according to a 2D histogram (in y)

Created on Mon 27Dec21: Version history:
----------------------------------------------
 1.0: 27Dec21: First implementation

 1.1: 22Jun21: Implemented differentiation of x and y to be able to generate phase space
               in the two different coordinates from two different histograms

@author: MarvinPfaff
"""

#--------  Module dependencies
import os
import ctypes
import ROOT

class RandomGenerator(object):

    __instance = None

#--------  "Built-in methods":
    def __new__(cls, rootfilename, histname, histname2Dx, histname2Dy):
        if cls.__instance is None:
            print('RandomGenerator.__new__: creating the RandomGenerator object')
            print('-------------------')
            cls.__instance = super(RandomGenerator, cls).__new__(cls)

            cls._rootfilename = rootfilename
            rootFile = ROOT.TFile(cls._rootfilename, 'READ', 'ROOT file with Histograms')
            cls._hist = rootFile.Get(histname)
            cls._hist.SetDirectory(0)
            cls._hist2Dx = rootFile.Get(histname2Dx)
            cls._hist2Dx.SetDirectory(0)
            cls._hist2Dy = rootFile.Get(histname2Dy)
            cls._hist2Dy.SetDirectory(0)
            rootFile.Close()

            # Summarise initialisation
            cls.print(cls)

        return cls.__instance

    def __repr__(self):
        return "RandomGenerator()"

    def __str__(self):
        self.__repr__()
        self.print()

    #--------  Module methods
    def getRandom(self):
        #rootFile = ROOT.TFile(self._rootfilename, 'READ', 'ROOT file with Histograms')
        #hist = rootFile.Get(histname)
        x = self._hist.GetRandom()
        return x

    def getRandom2Dx(self):
        x = ctypes.c_double(0.0)
        y = ctypes.c_double(0.0)
        #rootFile = ROOT.TFile(self._rootfilename, 'READ', 'ROOT file with Histograms')
        #hist2D = rootFile.Get(histname)
        self._hist2Dx.GetRandom2(x,y)
        return x.value, y.value

    def getRandom2Dy(self):
        x = ctypes.c_double(0.0)
        y = ctypes.c_double(0.0)
        #rootFile = ROOT.TFile(self._rootfilename, 'READ', 'ROOT file with Histograms')
        #hist2D = rootFile.Get(histname)
        self._hist2Dy.GetRandom2(x,y)
        return x.value, y.value

#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)

    def CdVrsn(self):
        return 1.1

    def getRootFilename(self):
        return self._rootfilename

#--------  Utilities:
    def print(self):
        print("    RandomGenerator.print: version:", self.CdVrsn(self))
        print("    ROOT filename for histogram input  :", self._rootfilename)
