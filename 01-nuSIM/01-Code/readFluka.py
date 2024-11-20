#       readFluka.py
#
#   Version 1.0                 Paul Kyberd
#
#   Read in events from the fluka file and pass the values for suitable events back to nuSIM
#
#	Instantiate with readFluka(meanP). meanP is the central value of the pion meomenta required
#	Particle with a PDFId of 211 (pi+) and a momentum meanP*0.9 < p < meanP*1.1 wil be returned
#   in value - with value[0]=x, value[1]=y, value[2]=px, value[3]=py, value[4]=pz. if there are
#	no further values an EOFError is raised (thrown)
#
#   Calling sequence
#       rf = readFluka.readFluka(meanP)
#   Can get the total number of entries in the file of all types using
#       nEntries = rf._entries
#	Read from the file into an array with all the data
#       rf.readEvent()
#	Total number of events in the array
#		nEve= numberOfEvents(self):
#   Then a suitable loop might be
#       for i in range(nEve):
#           values = rf.getNextEvent()
#
#	Error: if the getNextEvent() called passed the end of the array - the method exits with a message
#		exit("trying to read too many events")
#
#	Version history:
# ----------------------------------------------
#   1.0: 017Sept24: Initial version
#

import ROOT
import numpy as np

class readFluka:

	def __init__(self, meanP):
		self._curEvent = 0
		self._flukaFileName = "31-Target/Fluka100GeVHPos.root"
		self._flukaRootFile = ROOT.TFile(self._flukaFileName, 'read')
		self._ntuple = "Data"
		self._data = self._flukaRootFile.Get(self._ntuple)
		self._entries = self._data.GetEntries()
		self._meanP = meanP
		self._pLow = meanP*0.9
		self._pHigh = meanP*1.10
		print ("in readData data is ", self._data)
		self._values=[]
		self._eventsFound=0


	def readEvents(self):

		moreEvents = True
		self._eventsFound =0
		print ("entries is ", self._entries)
		while moreEvents:
			found = False
			while not found:
				self._data.GetEntry(self._curEvent)
				x = getattr(self._data, 'x')
				y = getattr(self._data, 'y')
				p = getattr(self._data, 'p')
				px = getattr(self._data, 'px')
				py = getattr(self._data, 'py')
				pz = getattr(self._data, 'pz')
				PDGId = getattr(self._data, 'PDGId')
				p = np.sqrt(px**2+py**2+pz**2)
				self._curEvent = self._curEvent + 1
				if ((PDGId == 211) and (p>=self._pLow) and (p<=self._pHigh)):
					found = True
					self._values.append([x,y,px,py,pz])
					self._eventsFound = self._eventsFound + 1
				if self._curEvent + 1 > self._entries:
					found = True
					self._flukaRootFile.Close()
					moreEvents = False

			if (self._curEvent%10000 ==0):
				print ("events found is ", self._eventsFound)
#				print ("curEvent is ", self._curEvent)

		print ("eventsFound is ", self._eventsFound)
		self._curEvent = 0
		return

	def numberOfEvents(self):
		return self._eventsFound

	def getNextEvent(self):
		if self._curEvent + 1 > self._eventsFound:
			exit("trying to read too many events")

		event = self._values[self._curEvent]
		self._curEvent = self._curEvent + 1
		return event


if __name__ == "__main__" : 

	print(" Does nothing when run as main")

