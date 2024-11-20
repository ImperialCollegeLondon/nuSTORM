"""
Storing and outputting record of a run
========================================

@date: 19 Nov 2024
@author: Paul Kyberd
@version: 1.0

@description: store and output summary information about a run

"""
class runRecord:

	dbg = True

	def __init__(self):
		if self.dbg: print("instantiating runRecord")

	def __str__(self):
		return f'runRecord: run {self._runNumber}'

	def setRunNumber(self, runNumber):
		self._runNumber = runNumber

	def openTextOutput(self, outputTextFileName):
		self._outputTextFile = open(outputTextFileName, "w")

	def writeTextFile(self):
		self._outputTextFile.write(f'run: {self._runNumber}')

	def closeTextOutput(self):
		self._outputTextFile.close()


if __name__ == "__main__" :
	rr = runRecord()
	rr.setRunNumber(34)
	rr.openTextOutput("outputTest.txt")
	rr.writeTextFile()
	rr.closeTextOutput()
	print(rr)