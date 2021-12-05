#   testUnit.py                                         Version 2.0
#       29 October 2021
#
#     Paul Kyberd                                       21 March 2011
#
#  Version 2.0                                          29 October 2021
#   Update to python 3
#
#  Version 1.0                                          21 March 2011
#   Some code to encapsulate tests
#

class testUnit():
    
    def __init__(self):
        self.testCount = 0
        self.failCount = 0
        self.Version = 1.0
    
    def announce(self,prog):
        print ("          Unit test harness: Version " + str(self.Version) + " --- testing " + str(prog))
        print ("  ")
  
    
    def fTry(self, routine):
        self.testCount = self.testCount + 1
#        print "fTry"
#        print routine[0]
#        print routine[1]
#        print routine[2]
        if routine[1] == routine[2]:
            print (routine[0] + ".... OK")
        else:
            print (routine[0] + ".... Failed")
            self.failCount = self.failCount + 1
            
#        print "end fTry"

    def summary(self):
        print (")Number of tests: " + str(self.testCount)+ "   Tests fails: " + str(self.failCount))

    def version(self):
        print ("Unit test harness: Version " + str(self.Version))
