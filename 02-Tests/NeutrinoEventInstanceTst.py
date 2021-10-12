#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "NeutrinoEventInstance" class
=============================================

  Assumes that nuSim code is in path.

  Script starts by testing built in methods.  Then a soak test with a
  large number of decays is executed.  Finally a set of reference plots
  are generated.

 03Apr21: Fix relativistic treatment of muon lifetime

"""

import os
import MuonConst as muCnst
import nuSTORMPrdStrght as nuPrdStrt
import NeutrinoEventInstance as nuEvtInst
import PionEventInstance as piEvtInst
import numpy as np
import matplotlib.pyplot as plt
import Simulation as Simu
import math as mth

mc = muCnst.MuonConst()

piEvt= piEvtInst.PionEventInstance(8)
tpi=piEvt.gettpi()
piTraceSpaceCoord=piEvt.getTraceSpaceCoord()
mu4mmtm=piEvt.getmu4mmtm()
mucostheta=piEvt.getcostheta()

##! Start:
print("========  NeutrinoEventInstance: tests start  ========")

debug=0

##! Create instance, test built-in methods:
NeutrinoEventInstanceTest = 1
print()
print("NeutrinoEventInstanceTest:", NeutrinoEventInstanceTest, \
      " Test built-in methods.")
nuSIMPATH = os.getenv('nuSIMPATH')
filename  = os.path.join(nuSIMPATH, \
                         '11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
                      
nuEI = nuEvtInst.NeutrinoEventInstance(tpi, piTraceSpaceCoord, mu4mmtm, mucostheta, filename)
print("    __str__: \n", nuEI)
print("    --repr__: \n", repr(nuEI))

##! Test get/set methods:
NeutrinoEventInstanceTest = 2
print()
print("NeutrinoEventInstanceTest:", NeutrinoEventInstanceTest, \
      "Test get/set methods.")
print("    Muon momentum (GeV):", nuEI.getpmu())

##! Test methods by which neutrino-creation event is generated:
NeutrinoEventInstanceTest = 3
print()
print("NeutrinoEventInstanceTest:", NeutrinoEventInstanceTest, \
      "Test methods by which neutrino-creation event is generated.")
nuStrt = nuPrdStrt.nuSTORMPrdStrght(filename)

x = nuEI.CreateNeutrinos(tpi,  piTraceSpaceCoord, mu4mmtm, mucostheta, nuStrt)
print("    Neutrino event: trace-space coordinates of muon at decay:")
print("    ----> P_e:", x[4])
print("    ----> P_nue:", x[5])
print("    ----> P_numu:", x[6])
del nuEI

##! Soak test:
NeutrinoEventInstanceTest = 4
print()
print("NeutrinoEventInstanceTest:", NeutrinoEventInstanceTest, \
      "soak test.")
#Pmu = 5.
#print("    ----> Muon momentum:", Pmu)
nuEI = []
piEI = []

for i in range(500):

    piEvt= piEvtInst.PionEventInstance(8)
    piEI.append(piEvt)
    tpi=piEvt.gettpi()
    piTraceSpaceCoord=piEvt.getTraceSpaceCoord()
    mu4mmtm=piEvt.getmu4mmtm()
    mucostheta=piEvt.getcostheta()
    nuEI.append(nuEvtInst.NeutrinoEventInstance(tpi,  piTraceSpaceCoord, mu4mmtm, mucostheta, filename))
for i in range(5):
    print("    piEI[",i, "]: \n", piEI[i])
    print("    nuEI[",i, "]: \n", nuEI[i])

##! Plot result of soak test:
NeutrinoEventInstanceTest = 5
print()
print("NeutrinoEventInstanceTest:", NeutrinoEventInstanceTest, \
      "plots from soak test.")

s       = np.array([])
x       = np.array([])
z       = np.array([])
Emu     = np.array([])
tmu     = np.array([])
xpi     = np.array([])
ypi     = np.array([])
pxpi    = np.array([])
pypi    = np.array([])

Ee      = np.array([])
Enue    = np.array([])
Enumu   = np.array([])
PBl     = []
Pe      = np.array([])     #.. Pe  stored in an np.array
Pel     = []               #.. Pe stored in a list for plotting
Pnue    = np.array([])
Pnuel   = []
Pnumu   = np.array([])
Pnumul  = []
tane    = np.array([])
tannue  = np.array([])
tannumu = np.array([])
Rnue    = np.array([])
Rnumu   = np.array([])
CutEnue = np.array([])
CutEnumu= np.array([])

Cose             = np.array([])  # Cosine of angles wrt beam direction
Cosnue           = np.array([])
Cosnumu          = np.array([])
PrdStrghtcose    = np.array([])
PrdStrghtcosnue  = np.array([])
PrdStrghtcosnumu = np.array([])
Arc1cose         = np.array([])
Arc1cosnue       = np.array([])
Arc1cosnumu      = np.array([])
RetStrghtcose    = np.array([])
RetStrghtcosnue  = np.array([])
RetStrghtcosnumu = np.array([])
Arc2cose         = np.array([])
Arc2cosnue       = np.array([])
Arc2cosnumu      = np.array([])

PrdStrghtLngth = nuStrt.ProdStrghtLen()
Circumference  = nuStrt.Circumference()
ArcLen         = nuStrt.ArcLen()
ArcRad         = ArcLen/mth.pi

for (piEvt,nuEvt) in zip(piEI,nuEI):
  #if(nuEvt.getAbsorbed()==False):
   # Case 1: for all decays except for the absorbed ones
  if(nuEvt.getAbsorbed()==False) and (nuEvt.getTraceSpaceCoord()[0]>PrdStrghtLngth):
  # Case 2: only for accepted muons, Not taking muon decays in the production straight before entering the ring into account
  # Case 3: Hash both if statements to get zero absorption
    xpi    = np.append(xpi,    piEvt.getTraceSpaceCoord()[1])
    ypi    = np.append(ypi,    piEvt.getTraceSpaceCoord()[2])
    pxpi   = np.append(pxpi,    piEvt.getmu4mmtm()[1][0]/piEvt.getmu4mmtm()[1][2])
    pypi   = np.append(pypi,    piEvt.getmu4mmtm()[1][1]/piEvt.getmu4mmtm()[1][2])
    
    tmu     = np.append(tmu,     nuEvt.gettmu())
    s     = np.append(s,     nuEvt.getTraceSpaceCoord()[0])
    where = nuEvt.getTraceSpaceCoord()[0]%Circumference
    
    xi    = nuEvt.getTraceSpaceCoord()[1]
    x     = np.append(x,     xi)
   
    zi    = nuEvt.getTraceSpaceCoord()[3]
    z     = np.append(z,     zi)
    
    Emu   = np.append(Emu,   piEvt.getmu4mmtm()[0])
    Ee    = np.append(Ee,    nuEvt.gete4mmtm()[0])
    Enue  = np.append(Enue,  nuEvt.getnue4mmtm()[0])
    Enumu = np.append(Enumu, nuEvt.getnumu4mmtm()[0])
    
    #3 momentum of the beam and decay products
    Pb     = nuEvt.getPb()
    PBl.append(Pb)
    Pb     = np.array(Pb)
    
    Pe3    = nuEvt.gete4mmtm()[1]
    Pel.append(Pe3)
    Pe    = np.append(Pe,    Pe3)
    
    Pnue3    = nuEvt.getnue4mmtm()[1]
    Pnuel.append(Pnue3)
    Pnue    = np.append(Pnue,    Pe3)
    
    Pnumu3  = nuEvt.getnumu4mmtm()[1]
    Pnumul.append(Pnumu3)
    Pnumu = np.append(Pnumu, Pnumu3)

    cose     = np.dot(Pb,Pe3)/(np.linalg.norm(Pb)*np.linalg.norm(Pe3))
    cosnue   = np.dot(Pb,Pnue3)/(np.linalg.norm(Pb)*np.linalg.norm(Pnue3))
    cosnumu  = np.dot(Pb,Pnumu3)/(np.linalg.norm(Pb)*np.linalg.norm(Pnumu3))
    
    Cose    = np.append(Cose, cose)
    Cosnue  = np.append(Cosnue, cosnue)
    Cosnumu = np.append(Cosnumu, cosnumu)
    
    if (PrdStrghtLngth>=where>=0):
        PrdStrghtcose = np.append(PrdStrghtcose,cose)
        PrdStrghtcosnue = np.append(PrdStrghtcosnue,cosnue)
        PrdStrghtcosnumu = np.append(PrdStrghtcosnumu,cosnumu)
    elif ((PrdStrghtLngth+ArcLen)>=where>PrdStrghtLngth):
        Arc1cose = np.append(Arc1cose,cose)
        Arc1cosnue = np.append(Arc1cosnue,cosnue)
        Arc1cosnumu = np.append(Arc1cosnumu,cosnumu)
    elif ((2*PrdStrghtLngth+ArcLen)>=where>(PrdStrghtLngth+ArcLen)):
        RetStrghtcose = np.append(RetStrghtcose,cose)
        RetStrghtcosnue = np.append(RetStrghtcosnue,cosnue)
        RetStrghtcosnumu = np.append(RetStrghtcosnumu,cosnumu)
    elif ((2*PrdStrghtLngth+2*ArcLen)>=where>(2*PrdStrghtLngth+ArcLen)):
        Arc2cose = np.append(Arc2cose,cose)
        Arc2cosnue = np.append(Arc2cosnue,cosnue)
        Arc2cosnumu = np.append( Arc2cosnumu,cosnumu)

    pt_e    = np.sqrt(nuEvt.gete4mmtm()[1][0]**2 \
                      + nuEvt.gete4mmtm()[1][1]**2)
    pt_nue  = np.sqrt(nuEvt.getnue4mmtm()[1][0]**2 \
                      + nuEvt.getnue4mmtm()[1][1]**2)
    pt_numu = np.sqrt(nuEvt.getnumu4mmtm()[1][0]**2 \
                      + nuEvt.getnumu4mmtm()[1][1]**2)

    tane    = np.append(tane,    (pt_e/nuEvt.gete4mmtm()[1][2]) )
    tannue  = np.append(tannue,  (pt_nue/nuEvt.getnue4mmtm()[1][2]) )
    tannumu = np.append(tannumu, (pt_numu/nuEvt.getnumu4mmtm()[1][2]) )
 
#-- Accepted Muon Trace Space Plots
plt.scatter(xpi,pxpi)
plt.xlabel('x')
plt.ylabel('xp')
plt.title('Accepted Muon Phase Space (X coordinates)')
plt.savefig('Scratch/NeutrinoEventInstanceTst_xp vs x.pdf')
#plt.show()
plt.close()

plt.scatter(ypi,pypi)
plt.xlabel('y')
plt.ylabel('yp')
plt.title('Accepted Muon Phase Space (Y coordinates)')
plt.savefig('Scratch/NeutrinoEventInstanceTst_yp vs y.pdf')
#plt.show()
plt.close()

plt.scatter(pxpi,Emu)
plt.xlabel('xp')
plt.ylabel('Muon Energy Emu (Gev)')
plt.title('Muon energy vs xp plot')
plt.savefig('Scratch/NeutrinoEventInstanceTst_Emu vs xp.pdf')
#plt.show()
plt.close()    

plt.scatter(pypi,Emu)
plt.xlabel('yp')
plt.ylabel('Muon Energy Emu (Gev)')
plt.title('Muon energy vs yp')
plt.savefig('Scratch/NeutrinoEventInstanceTst_Emu vs yp.pdf')
#plt.show()
plt.close()

#-- Lifetime distribution:
n, bins, patches = plt.hist(s, bins=50, color='y', log=True)
plt.xlabel('s (m)')
plt.ylabel('Entries')
plt.title('s distribution')
# add a 'best fit' line  Error: Pmu not defined because no central muon momentum is defined 
#beta  = Pmu/(mth.sqrt(Pmu**2 + (mc.mass()/1000.)**2))
#gamma = (mth.sqrt(Pmu**2 + (mc.mass()/1000.)**2))/(mc.mass()/1000.)
#l = 1./(gamma*mc.lifetime()*beta*mc.SoL())
#y = n[0]*np.exp(-l*bins)
#plt.plot(bins, y, '-', color='b')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot1.pdf')
#plt.show()
plt.close()

#-- Lifetime distribution:
n, bins, patches = plt.hist(tmu, bins=50, color='y', log=True)
plt.xlabel('tmu (s)')
plt.ylabel('Entries')
plt.title('Muon Decay Time Distribution')
# add a 'best fit' line  
l = 1/mc.lifetime()
y = n[0]*np.exp(-l*bins)
plt.plot(bins, y, '-', color='b')
plt.savefig('Scratch/NeutrinoEventInstanceTst_tmu.pdf')
plt.show()
plt.close()

#-- Energy distributions:

n, bins, patches = plt.hist(Emu, bins=50, color='y', range=(0,9))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('Muon energy distribution')
plt.savefig('Scratch/NeutrinoEventInstanceTst_Muon Energy Distribution.pdf')
plt.show()
plt.close()

n, bins, patches = plt.hist(Ee, bins=50, color='y', range=(0.,5.0))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('Electron energy distribution')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot2.pdf')
#plt.show()
plt.close()

n, bins, patches = plt.hist(Enue, bins=50, color='y', range=(0.,5.))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('Electron-neutrino energy distribution')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot3.pdf')
#plt.show()
plt.close()

n, bins, patches = plt.hist(Enumu, bins=50, color='y', range=(0.,5.))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('Muon-neutrino energy distribution')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot4.pdf')
#plt.show()
plt.close()

## Z vs X distribution
plt.scatter(z,x)
plt.axis([max(z)+20,min(z)-20,max(x)+20,min(x)-20])
plt.xlabel('Z')
plt.ylabel('X')
plt.title('Muon Decay Position')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot5.pdf')
#plt.show()
plt.close()

##beamdirectionplotting
fig, ax = plt.subplots()
U=np.array([])
V=np.array([])
i=1
for P in PBl:
    if(debug==1):
          print("Plotting: Beam momentum Pb",P)
          print("Number of instance:", i)
          i=i+1
    U=np.append(U,- P[2])
    V=np.append(V,- P[0])
    # X= nuSTORM -Z axis; # Y= nuSTORM -X axis
Q=ax.quiver(z,x, U, V)
ax.quiverkey(Q, X=0.5, Y=0.90, U=5, label='Beam Momentum', labelpos='E')
ax.invert_xaxis()
ax.invert_yaxis()
plt.xlabel('Z')
plt.ylabel('X')
plt.title('Beam Direction')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot18.pdf')
plt.show()
plt.close()


##Decay product momentum direction
fig, ax = plt.subplots()
U=np.array([])
V=np.array([])
i=1
for P in PBl:
    if(debug==1):
          print("Plotting: Beam momentum Pb",P)
          print("Number of instance:", i)
          i=i+1
    U=np.append(U,- P[2]) 
    V=np.append(V,- P[0])  
    # X= nuSTORM -Z axis; # Y= nuSTORM -X axis
Q=ax.quiver(z,x, U, V)  
ax.quiverkey(Q, X=0.4, Y=0.90, U=5, label='Beam Momentum', labelpos='E')
U=np.array([])
V=np.array([])
i=1
for P in Pel:
    
    if(debug==1):
          print("Plotting: electron momentum Pel",P)
          print("Number of instance:", i)
          i=i+1
    U=np.append(U,- P[2])
    V=np.append(V,- P[0])
Q=ax.quiver(z,x, U, V, color='r')
ax.quiverkey(Q, X=0.4, Y=0.85, U=2, label='Electron Momentum', labelpos='E')
ax.invert_xaxis()
ax.invert_yaxis()
plt.xlabel('Z')
plt.ylabel('X')
plt.title('Electron Momentum Vectors')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot23.pdf')
plt.close()
fig, ax = plt.subplots()
U=np.array([])
V=np.array([])
i=1
for P in PBl:
    if(debug==1):
          print("Plotting: Beam momentum Pb",P)
          print("Number of instance:", i)
          i=i+1
    U=np.append(U,- P[2])
    V=np.append(V,- P[0])
    # X= nuSTORM -Z axis; # Y= nuSTORM -X axis
Q=ax.quiver(z,x, U, V)
ax.quiverkey(Q, X=0.4, Y=0.90, U=5, label='Beam Momentum', labelpos='E')
U=np.array([])
V=np.array([])
i=1
for P in Pnuel:
    if(debug==1):
          print("Plotting: electron neutrino momentum Pnuel",P)
          print("Number of instance:", i)
          i=i+1
    U=np.append(U,- P[2])
    V=np.append(V,- P[0])
Q=ax.quiver(z,x, U, V, color='g')
ax.quiverkey(Q, X=0.4, Y=0.80, U=1.5, label='Electron Neutrino Momentum', labelpos='E')
ax.invert_xaxis()
ax.invert_yaxis()
plt.xlabel('Z')
plt.ylabel('X')
plt.title('Electron Neutrino Momentum Vectors')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot24.pdf')
plt.close()

fig, ax = plt.subplots()
U=np.array([])
V=np.array([])
i=1
for P in PBl:
    if(debug==1):
          print("Plotting: Beam momentum Pb",P)
          print("Number of instance:", i)
          i=i+1
    U=np.append(U,- P[2])
    V=np.append(V,- P[0])
    # X= nuSTORM -Z axis; # Y= nuSTORM -X axis
Q=ax.quiver(z,x, U, V)
ax.quiverkey(Q, X=0.4, Y=0.90, U=5, label='Beam Momentum', labelpos='E')

U=np.array([])
V=np.array([])
i=1

for P in Pnumul:
    if(debug==1):
          print("Plotting: muon neutrino momentum Pnumul",P)
          print("Number of instance:", i)
          i=i+1
    
    U=np.append(U,- P[2])
    V=np.append(V,- P[0])
Q=ax.quiver(z,x, U, V, color='b')
ax.quiverkey(Q, X=0.4, Y=0.75, U=2, label='Muon Neutrino Momentum', labelpos='E')
ax.invert_xaxis()
ax.invert_yaxis()
plt.xlabel('Z')
plt.ylabel('X')
plt.title('Muon Neutrino Momentum Vectors')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot25.pdf')
plt.close()


##Angular distribution of decay product momentum with respect to beam momentum
n, bins, patches = plt.hist(PrdStrghtcose, bins=500, color='y', range=(0.99,1))
plt.xlabel('PrdStrghtcose')
plt.ylabel('Frequency')
plt.title('Angular distribution of electron momentum \n with respect to beam direction in the production straight')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot6.pdf')
plt.close()

n, bins, patches = plt.hist(PrdStrghtcosnue, bins=500, color='y', range=(0.99,1))
plt.xlabel('PrdStrghtcosnue')
plt.ylabel('Frequency')
plt.title('Angular distribution of electron neutrino momentum \n with respect to beam direction in the production straight')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot7.pdf')
plt.close()

n, bins, patches = plt.hist(PrdStrghtcosnumu, bins=500, color='y', range=(0.99,1))
#print(PrdStrghtcosnumu)
plt.xlabel('PrdStrghtcosnumu')
plt.ylabel('Frequency')
plt.title('Angular distribution of muon neutrino momentum \n with respect to beam direction in the production straight')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot8.pdf')
plt.close()

n, bins, patches = plt.hist(Arc1cose, bins=500, color='y', range=(0.99,1))
plt.xlabel('Arc1cose')
plt.ylabel('Frequency')
plt.title('Angular distribution of electron momentum \n with respect to beam direction in the first arc')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot9.pdf')
plt.close()

n, bins, patches = plt.hist(Arc1cosnue, bins=500, color='y', range=(0.99,1))
plt.xlabel('Arc1cosnue')
plt.ylabel('Frequency')
plt.title('Angular distribution of electron neutrino momentum \n with respect to beam direction in the first arc')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot10.pdf')
plt.close()

n, bins, patches = plt.hist(Arc1cosnumu, bins=500, color='y', range=(0.99,1))
#print(Arc1cosnumu)
plt.xlabel('Arc1cosnumu')
plt.ylabel('Frequency')
plt.title('Angular distribution of muon neutrino momentum \n with respect to beam direction in the first arc')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot11.pdf')
plt.close()

n, bins, patches = plt.hist(RetStrghtcose, bins=500, color='y', range=(0.99,1))
plt.xlabel('RetStrghtcose')
plt.ylabel('Frequency')
plt.title('Angular distribution of electron momentum \n with respect to beam direction in the return straight')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot12.pdf')
plt.close()

n, bins, patches = plt.hist(RetStrghtcosnue, bins=500, color='y', range=(0.99,1))
plt.xlabel('RetStrghtcosnue')
plt.ylabel('Frequency')
plt.title('Angular distribution of electron neutrino momentum \n with respect to beam direction in the return straight')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot13.pdf')
plt.close()

n, bins, patches = plt.hist(RetStrghtcosnumu, bins=500, color='y', range=(0.99,1))
#print(RetStrghtcosnumu)
plt.xlabel('RetStrghtcosnumu')
plt.ylabel('Frequency')
plt.title('Angular distribution of electron neutrino momentum \n with respect to beam direction in the return straight')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot14.pdf')
plt.close()

n, bins, patches = plt.hist(Arc2cose, bins=500, color='y', range=(0.99,1))
plt.xlabel('Arc2cose')
plt.ylabel('Frequency')
plt.title('Angular distribution of electron momentum \n with respect to beam direction in the second arc')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot15.pdf')
plt.close()

n, bins, patches = plt.hist(Arc2cosnue, bins=500, color='y', range=(0.99,1))
plt.xlabel('Arc2cosnue')
plt.ylabel('Frequency')
plt.title('Angular distribution of electron neutrino momentum \n with respect to beam direction in the second arc')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot16.pdf')
plt.close()

#print(Arc2cosnumu)
n, bins, patches = plt.hist(Arc2cosnumu, bins=500, color='y', range=(0.99,1))
plt.xlabel('Arc2cosnumu')
plt.ylabel('Frequency')
plt.title('Angular distribution of muon neutrino momentum \n with respect to beam direction in the second arc')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot17.pdf')
plt.close()

n, bins, patches = plt.hist(Cose, bins=500, color='r', range=(0.99,1))
plt.xlabel('Cos(Angle between electron direction and beam direction)')
plt.ylabel('Frequency')
plt.title('Angular distribution of electron momentum \n with respect to beam direction')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot20.pdf')
plt.close()

n, bins, patches = plt.hist(Cosnue, bins=500, color='g', range=(0.99,1))
plt.xlabel('Cos(Angle between electron neutrino direction and beam direction)')
plt.ylabel('Frequency')
plt.title('Angular distribution of electron neutrino momentum \n with respect to beam direction')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot21.pdf')
plt.close()

n, bins, patches = plt.hist(Cosnumu, bins=500, color='b', range=(0.99,1))
plt.xlabel('Cos(Angle between muon neutrino direction and beam direction)')
plt.ylabel('Frequency')
plt.title('Angular distribution of muon neutrino momentum \n with respect to beam direction')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot22.pdf')
plt.show()
plt.close()



##! Complete:
print()
print("========  NeutrinoEventInstance: tests complete  ========")
