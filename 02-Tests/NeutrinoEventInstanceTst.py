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
import numpy as np
import matplotlib.pyplot as plt
import Simulation as Simu
import math as mth

mc = muCnst.MuonConst()

##! Start:
print("========  NeutrinoEventInstance: tests start  ========")

debug=0

##! Create instance, test built-in methods:
NeutrinoEventInstanceTest = 1
print()
print("NeutrinoEventInstanceTest:", NeutrinoEventInstanceTest, " Test built-in methods.")
nuSIMPATH = os.getenv('nuSIMPATH')
filename  = os.path.join(nuSIMPATH, '11-Parameters/nuSTORM-PrdStrght-Params-v1.0.csv')
nuEI = nuEvtInst.NeutrinoEventInstance(5., filename)
print("    __str__:", nuEI)
print("    --repr__", repr(nuEI))

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
x = nuEI.CreateNeutrinos(nuStrt)
print("    Neutrino event: trace-space coordinates of muon at decay, P_e, P_nue, P_numu:", x)
del nuEI

##! Soak test:
NeutrinoEventInstanceTest = 4
print()
print("NeutrinoEventInstanceTest:", NeutrinoEventInstanceTest, \
      "soak test.")
Pmu = 5.
print("    ----> Muon momentum:", Pmu)
nuEI = []
for i in range(10000):
    if(debug==1):
      print("Number of instance:", i+1)
    nuEI.append(nuEvtInst.NeutrinoEventInstance(Pmu, filename))
for i in range(5):
    print("    nuEI[i]:", nuEI[i])

##! Plot result of soak test:
NeutrinoEventInstanceTest = 5
print()
print("NeutrinoEventInstanceTest:", NeutrinoEventInstanceTest, \
      "plots from soak test.")

s       = np.array([])
x       = np.array([])
z       = np.array([])
Ee      = np.array([])
Enue    = np.array([])
Enumu   = np.array([])
PBl     = []
Pe      = np.array([])
Pel     = []
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

if (debug==1):
    print("Soak Test:PrdStrghtLngth, Circumference, ArcLen, ArcRad ",PrdStrghtLngth, Circumference, ArcLen, ArcRad)
    i=1
for nuEvt in nuEI:
    if (debug==1):
         print("Number of instance:", i)
         i=i+1
    s     = np.append(s,     nuEvt.getTraceSpaceCoord()[0])
    where = nuEvt.getTraceSpaceCoord()[0]%Circumference
    
    xi    = nuEvt.getTraceSpaceCoord()[1]
    x     = np.append(x,     xi)
   
    zi    = nuEvt.getTraceSpaceCoord()[3]
    z     = np.append(z,     zi)
    
    if(debug==1):
        print("Soak Test: s, where, xi, zi",nuEvt.getTraceSpaceCoord()[0], where, xi, zi)
        
   

 

    Ee    = np.append(Ee,    nuEvt.gete4mmtm()[0])
    Enue  = np.append(Enue,  nuEvt.getnue4mmtm()[0])
    Enumu = np.append(Enumu, nuEvt.getnumu4mmtm()[0])
    
    #3 momentum of the beam and decay products
    Pb     = nuEvt.getPb()
    if(debug==1):
      print("SoakTest: Pb", Pb)
    PBl.append(Pb)
    Pb     = np.array(Pb)
    
 
    Pe3    = nuEvt.gete4mmtm()[1]
    if(debug==1):
      print("SoakTest: Pe3", Pe3)
    Pel.append(Pe3)
    Pe    = np.append(Pe,    Pe3)
    
    
    Pnue3    = nuEvt.getnue4mmtm()[1]
    if(debug==1):
      print("SoakTest: Pnue3", Pnue3)
    Pnuel.append(Pnue3)
    Pnue    = np.append(Pnue,    Pe3)
    
    
    Pnumu3  = nuEvt.getnumu4mmtm()[1]
    if(debug==1):
      print("SoakTest: Pnumu3", Pnumu3)
    Pnumul.append(Pnumu3)
    Pnumu = np.append(Pnumu, Pnumu3)
   

    #cosine of angle between the beam direction and resultant product direction
    cose     = np.dot(Pb,Pe3)/(np.linalg.norm(Pb)*np.linalg.norm(Pe3))
    cosnue   = np.dot(Pb,Pnue3)/(np.linalg.norm(Pb)*np.linalg.norm(Pnue3))
    cosnumu  = np.dot(Pb,Pnumu3)/(np.linalg.norm(Pb)*np.linalg.norm(Pnumu3))
    

    #if(PrdStrghtLngth>=zi>=0) and (1>=xi>=-1):
    if (PrdStrghtLngth>=where>=0):
        if(debug==1):
          print("Soak Test: In the production straight")
        PrdStrghtcose = np.append(PrdStrghtcose,cose)
        PrdStrghtcosnue = np.append(PrdStrghtcosnue,cosnue)
        PrdStrghtcosnumu = np.append(PrdStrghtcosnumu,cosnumu)
    if ((PrdStrghtLngth+ArcLen)>=where>PrdStrghtLngth):
    #if(zi>PrdStrghtLngth):
        if(debug==1):
          print("Soak Test: In the Arc 1")
        Arc1cose = np.append(Arc1cose,cose)
        Arc1cosnue = np.append(Arc1cosnue,cosnue)
        Arc1cosnumu = np.append(Arc1cosnumu,cosnumu)
    if ((2*PrdStrghtLngth+ArcLen)>=where>(PrdStrghtLngth+ArcLen)):
        if(debug==1):
          print("Soak Test: In the return straight")
        RetStrghtcose = np.append(RetStrghtcose,cose)
        RetStrghtcosnue = np.append(RetStrghtcosnue,cosnue)
        RetStrghtcosnumu = np.append(RetStrghtcosnumu,cosnumu)
   
    if ((2*PrdStrghtLngth+2*ArcLen)>=where>(2*PrdStrghtLngth+ArcLen)):
        if(debug==1):
          print("Soak Test: In the return straight")
        Arc2cose = np.append(Arc2cose,cose)
        Arc2cosnue = np.append(Arc2cosnue,cosnue)
        Arc2cosnumu = np.append( Arc2cosnumu,cosnumu)
  
    pt_e    = np.sqrt(nuEvt.gete4mmtm()[1][0]**2 + nuEvt.gete4mmtm()[1][1]**2)
    pt_nue  = np.sqrt(nuEvt.getnue4mmtm()[1][0]**2 + nuEvt.getnue4mmtm()[1][1]**2)
    pt_numu = np.sqrt(nuEvt.getnumu4mmtm()[1][0]**2 + nuEvt.getnumu4mmtm()[1][1]**2)

    tane    = np.append(tane,    (pt_e/nuEvt.gete4mmtm()[1][2]) )
    tannue  = np.append(tannue,  (pt_nue/nuEvt.getnue4mmtm()[1][2]) )
    tannumu = np.append(tannumu, (pt_numu/nuEvt.getnumu4mmtm()[1][2]) )


#-- Lifetime distribution:
n, bins, patches = plt.hist(s, bins=50, color='y', log=True)
plt.xlabel('s (m)')
plt.ylabel('Entries')
plt.title('s distribution')
# add a 'best fit' line
Emu   = mth.sqrt(Pmu**2 + (mc.mass()/1000.)**2)
beta  = Pmu/Emu
gamma = Emu/(mc.mass()/1000.)
l = 1./(gamma*mc.lifetime()*beta*mc.SoL())
y = n[0]*np.exp(-l*bins)
plt.plot(bins, y, '-', color='b')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot1.pdf')
plt.close()

#-- Energy distributions:
n, bins, patches = plt.hist(Ee, bins=50, color='y', range=(0.,5.0))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('Electron energy distribution')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot2.pdf')
plt.close()

n, bins, patches = plt.hist(Enue, bins=50, color='y', range=(0.,5.))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('Electron-neutrino energy distribution')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot3.pdf')
plt.close()

n, bins, patches = plt.hist(Enumu, bins=50, color='y', range=(0.,5.))
plt.xlabel('Energy (GeV)')
plt.ylabel('Frequency')
plt.title('Muon-neutrino energy distribution')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot4.pdf')
plt.close()

## Z vs X distribution
plt.scatter(z,x)
plt.axis([max(z)+20,min(z)-20,max(x)+20,min(x)-20])
plt.xlabel('Z')
plt.ylabel('X')
plt.title('Muon Decay Position')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot5.pdf')
plt.close()

debug=0
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
plt.title('Momentum of Decay Products')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot19.pdf')
plt.close()


##Angular distribution of decay product momentum with respect to beam momentum
n, bins, patches = plt.hist(PrdStrghtcose, bins=50, color='y', range=(min(PrdStrghtcose),max(PrdStrghtcose)))
plt.xlabel('PrdStrghtcose')
plt.ylabel('Frequency')
plt.title('Angular distrubution of electron momentum \n with respect to beam direction in the production straight')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot6.pdf')
plt.close()

n, bins, patches = plt.hist(PrdStrghtcosnue, bins=50, color='y', range=(min(PrdStrghtcosnue),max(PrdStrghtcosnue)))
plt.xlabel('PrdStrghtcosnue')
plt.ylabel('Frequency')
plt.title('Angular distrubution of electron neutrino momentum \n with respect to beam direction in the production straight')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot7.pdf')
plt.close()

n, bins, patches = plt.hist(PrdStrghtcosnumu, bins=50, color='y', range=(min(PrdStrghtcosnumu),max(PrdStrghtcosnumu)))
#print(PrdStrghtcosnumu)
plt.xlabel('PrdStrghtcosnumu')
plt.ylabel('Frequency')
plt.title('Angular distrubution of muon neutrino momentum \n with respect to beam direction in the production straight')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot8.pdf')
plt.close()

n, bins, patches = plt.hist(Arc1cose, bins=50, color='y', range=(min(Arc1cose),max(Arc1cose)))
plt.xlabel('Arc1cose')
plt.ylabel('Frequency')
plt.title('Angular distrubution of electron momentum \n with respect to beam direction in the first arc')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot9.pdf')
plt.close()

n, bins, patches = plt.hist(Arc1cosnue, bins=50, color='y', range=(min(Arc1cosnue),max(Arc1cosnue)))
plt.xlabel('Arc1cosnue')
plt.ylabel('Frequency')
plt.title('Angular distrubution of electron neutrino momentum \n with respect to beam direction in the first arc')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot10.pdf')
plt.close()

n, bins, patches = plt.hist(Arc1cosnumu, bins=50, color='y', range=(min(Arc1cosnumu),max(Arc1cosnumu)))
#print(Arc1cosnumu)
plt.xlabel('Arc1cosnumu')
plt.ylabel('Frequency')
plt.title('Angular distrubution of muon neutrino momentum \n with respect to beam direction in the first arc')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot11.pdf')
plt.close()

n, bins, patches = plt.hist(RetStrghtcose, bins=50, color='y', range=(min(RetStrghtcose),max(RetStrghtcose)))
plt.xlabel('RetStrghtcose')
plt.ylabel('Frequency')
plt.title('Angular distrubution of electron momentum \n with respect to beam direction in the production straight')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot12.pdf')
plt.close()

n, bins, patches = plt.hist(RetStrghtcosnue, bins=50, color='y', range=(min(RetStrghtcosnue),max(RetStrghtcosnue)))
plt.xlabel('RetStrghtcosnue')
plt.ylabel('Frequency')
plt.title('Angular distrubution of electron neutrino momentum \n with respect to beam direction in the return straight')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot13.pdf')
plt.close()

n, bins, patches = plt.hist(RetStrghtcosnumu, bins=50, color='y', range=(min(RetStrghtcosnumu),max(RetStrghtcosnumu)))
#print(RetStrghtcosnumu)
plt.xlabel('RetStrghtcosnumu')
plt.ylabel('Frequency')
plt.title('Angular distrubution of electron neutrino momentum \n with respect to beam direction in the return straight')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot14.pdf')
plt.close()

n, bins, patches = plt.hist(Arc2cose, bins=50, color='y', range=(min(Arc2cose),max(Arc2cose)))
plt.xlabel('Arc2cose')
plt.ylabel('Frequency')
plt.title('Angular distrubution of electron momentum \n with respect to beam direction in the second arc')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot15.pdf')
plt.close()

n, bins, patches = plt.hist(Arc2cosnue, bins=50, color='y', range=(min(Arc2cosnue),max(Arc2cosnue)))
plt.xlabel('Arc2cosnue')
plt.ylabel('Frequency')
plt.title('Angular distrubution of electron neutrino momentum \n with respect to beam direction in the second arc')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot16.pdf')
plt.close()

n, bins, patches = plt.hist(Arc2cosnumu, bins=50, color='y', range=(min(Arc2cosnumu),max(Arc2cosnumu)))
#print(Arc2cosnumu)
plt.xlabel('Arc2cosnumu')
plt.ylabel('Frequency')
plt.title('Angular distrubution of muon neutrino momentum \n with respect to beam direction in the second arc')
plt.savefig('Scratch/NeutrinoEventInstanceTst_plot17.pdf')
plt.close()
##! Complete:
print()
print("========  NeutrinoEventInstance: tests complete  ========")
