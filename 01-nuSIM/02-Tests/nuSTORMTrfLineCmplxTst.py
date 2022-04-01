#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "nuSTORMTrfLineCmplx" class
========================================

  Assumes python path includes nuSim code.

  Script starts by testing built in methods.  Then a soak test with a set
  of reference plots.

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from iminuit import Minuit
from iminuit.cost import LeastSquares
from Simulation import *
import nuSTORMTrfLineCmplx as nuTrf
import PionEventInstance as PionEventInstance
import eventHistory
import PionConst as PionConst
import MuonConst as MuonConst

piCnst = PionConst.PionConst()
muCnst = MuonConst.MuonConst()

n_tst = 200000
t_bunch = 2.
t_spacing = 5.
t_extraction = 10.5

##! Start:
print("========  nuSTORMTrfLineCmplx: tests start  ========")

##! Test singleton class feature:
nuSTORMTrfLineCmplxTest = 1
print()
print("nuSTORMTrfLineCmplxTest:", nuSTORMTrfLineCmplxTest, " check if class is a singleton.")
nuSIMPATH = os.getenv('nuSIMPATH')
filename  = os.path.join(nuSIMPATH,r'11-Parameters/nuSTORM-TrfLineCmplx-Params-v1.0.csv')
#filename = r'/home/marvin/Documents/masters_thesis/nuSTORM/11-Parameters/nuSTORM-TrfLine-Params-v1.0.csv'
nuTrfLineCmplx  = nuTrf.nuSTORMTrfLineCmplx(filename)
nuTrfLineCmplx1 = nuTrf.nuSTORMTrfLineCmplx(filename)
print("    nuTrfLineCmplx singleton test:", id(nuTrfLineCmplx), id(nuTrfLineCmplx1), id(nuTrfLineCmplx)-id(nuTrfLineCmplx1))
if nuTrfLineCmplx != nuTrfLineCmplx1:
    raise Exception("nuSTORMTrfLineCmplx is not a singleton class!")

##! Check built-in methods:
nuSTORMTrfLineCmplxTest = 2
print()
print("nuSTORMTrfLineCmplxTest:", nuSTORMTrfLineCmplxTest, " check built-in methods.")
print("    __repr__:")
print(nuTrfLineCmplx)

##! Check get methods:
nuSTORMTrfLineCmplxTest = 3
print()
print("nuSTORMTrfLineCmplxTest:", nuSTORMTrfLineCmplxTest, " check get methods.")
print("----> printParams() method; reports parameters loaded")
nuTrfLineCmplx.printParams()

##! Check momentum, z, and transverse distributions:
nuSTORMTrfLineCmplxTest = 4
print()
print("nuSTORMTrfLineCmplxTest:", nuSTORMTrfLineCmplxTest, " momentum distribution.")
n_tst4 = n_tst
npi = 0
ppi = np.array([])
Epi = np.array([])
tpi = np.array([])
zpi = np.array([])
xpi = np.array([])
ypi = np.array([])
xppi = np.array([])
yppi = np.array([])
t_dcy_pi = np.array([])
s_dcy_pi = np.array([])
t_prod_pi = np.array([])
print()
for i in range(n_tst4):
  npi += 1
  ppi = np.append(ppi, nuTrfLineCmplx.GeneratePiMmtm(5.))
  Epi = np.append(Epi, nuTrfLineCmplx.CalculateE(ppi[npi-1]))
  tpi = np.append(tpi, nuTrfLineCmplx.GenerateTime())
  spi = tpi[npi-1]*nuTrfLineCmplx.Calculatev(ppi[npi-1])
  zpi = np.append(zpi, nuTrfLineCmplx.TrfLineCmplxLen()+nuTrfLineCmplx.Calculatez(spi))
  x, y, xp, yp = nuTrfLineCmplx.GenerateTrans(zpi[npi-1])
  xpi  = np.append(xpi, x)
  ypi  = np.append(ypi, y)
  xppi = np.append(xppi, xp)
  yppi = np.append(yppi, yp)
  gamma = Epi[i]/(piCnst.mass()/1000.)
  vpi = nuTrfLineCmplx.Calculatev(ppi[npi-1])
  t_dcy_pi = np.append(t_dcy_pi, nuTrfLineCmplx.GenerateDcyTime(gamma))
  DcyPt = nuTrfLineCmplx.CalculateDcyPt(s=0.,s_final=nuTrfLineCmplx.TrfLineCmplxLen(),t_i=tpi[npi-1],v=vpi,gamma=gamma)
  if DcyPt > -100.:
      s_dcy_pi = np.append(s_dcy_pi, DcyPt)
  t_prod_pi = np.append(t_prod_pi, nuTrfLineCmplx.Calculatet(s=0.,s_final=nuTrfLineCmplx.TrfLineCmplxLen(),t_i=tpi[npi-1],v=vpi))
  if npi < 11:
      print(npi, ppi[npi-1], tpi[npi-1], zpi[npi-1], x, y, xp, yp)

'''
##! Check generating pion distributions from scratch
nuSTORMTrfLineCmplxTest = 5
print()
print("nuSTORMTrfLineCmplxTest:", nuSTORMTrfLineCmplxTest, " pion distributions from scratch.")
p0 = 5.
runNum = 101
eventNum = 0
n_tst5 = n_tst
weight = 1/n_tst5

dir = os.path.join(nuSIMPATH,r'Scratch/nuSTORMTrfLineCmplx')
dirExist = os.path.isdir(dir)
if dirExist == False:
    os.mkdir(dir)

eH = eventHistory.eventHistory()
eH.outFile("Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst5.root")
eH.rootStructure()

#pion at target
x_pi_t = np.array([])
y_pi_t = np.array([])
E_pi_t = np.array([])
xp_pi_t = np.array([])
yp_pi_t = np.array([])
t_pi_t = np.array([])

#pion at production straight entrance
t_pi_prod = np.array([])

#pion at pion decay
s_pi_dcy = np.array([])
t_pi_dcy = np.array([])

#muon from pion decay
x_mu = np.array([])
y_mu = np.array([])
E_mu = np.array([])
xp_mu = np.array([])
yp_mu = np.array([])

#neutrino from pion decay
x_nu = np.array([])
y_nu = np.array([])
E_nu = np.array([])
xp_nu = np.array([])
yp_nu = np.array([])

for eventNum in range(n_tst5):

    #this is just needed for testing in order to reset the event history at the beginning of each loop
    #(otherwise locations that have not been filled this time but before will retain the earlier values)
    #eH = eventHistory.eventHistory()

    eH = nuTrfLineCmplx.GeneratePion(eventHist=eH,runNum=runNum,eventNum=eventNum,weight=weight,p0=p0)

    pi_target = eH.findParticle('target')
    x_pi_t = np.append(x_pi_t, pi_target.traceSpace().x())
    y_pi_t = np.append(y_pi_t, pi_target.traceSpace().y())
    E_pi_t = np.append(E_pi_t, pi_target.p()[0])
    xp_pi_t = np.append(xp_pi_t, pi_target.traceSpace().xp())
    yp_pi_t = np.append(yp_pi_t, pi_target.traceSpace().yp())
    t_pi_t = np.append(t_pi_t, pi_target.t())

    pi_prod = eH.findParticle('productionStraight')
    pi_dcy = eH.findParticle('pionDecay')
    mu = eH.findParticle('muonProduction')
    nu = eH.findParticle('piFlashNu')

    eH.fill()

    if pi_prod.run() != -1:
        t_pi_prod = np.append(t_pi_prod, pi_prod.t())

        if pi_dcy.run() != -1 or mu.run() != -1 or nu.run() != -1:
            raise Exception("ERROR: Pion at production straight entrance generated even though pion decayed!")

    elif pi_dcy.run() != -1:
        s_pi_dcy = np.append(s_pi_dcy, pi_dcy.traceSpace().s())
        t_pi_dcy = np.append(t_pi_dcy, pi_dcy.t())

        x_mu = np.append(x_mu, mu.traceSpace().x())
        y_mu = np.append(y_mu, mu.traceSpace().y())
        E_mu = np.append(E_mu, mu.p()[0])
        xp_mu = np.append(xp_mu, mu.traceSpace().xp())
        yp_mu = np.append(yp_mu, mu.traceSpace().yp())

        x_nu = np.append(x_nu, nu.traceSpace().x())
        y_nu = np.append(y_nu, nu.traceSpace().y())
        E_nu = np.append(E_nu, nu.p()[0])
        xp_nu = np.append(xp_nu, nu.traceSpace().xp())
        yp_nu = np.append(yp_nu, nu.traceSpace().yp())

        if pi_prod.run() != -1:
            raise Exception("ERROR: Pion at production straight entrance generated even though pion decayed!")

        if mu.run() == -1 or nu.run() == -1:
            raise Exception("ERROR: Pion decayed but no daughter particle(s) generated.")

eH.write()
eH.outFileClose()

##! Check generating pion distributions from parameter input
nuSTORMTrfLineCmplxTest = 6
print()
print("nuSTORMTrfLineCmplxTest:", nuSTORMTrfLineCmplxTest, " pion distributions from parameter input.")
runNum = 102
eventNum = 0
n_tst6 = n_tst4
weight = 1/n_tst6

n_dcy = 0

E_pi_trgt_tst6 = np.array([])

t_pi_prod_tst6 = np.array([])

t_pi_dcy_tst6 = np.array([])
z_pi_dcy_tst6 = np.array([])

E_mu_tst6 = np.array([])
x_mu_tst6 = np.array([])
y_mu_tst6 = np.array([])
xp_mu_tst6 = np.array([])
yp_mu_tst6 = np.array([])

E_mu_rf_tst6 = np.array([])
xp_mu_rf_tst6 = np.array([])
yp_mu_rf_tst6 = np.array([])

cosTheta = np.array([])
phi = np.array([])

E_nu_tst6 = np.array([])
x_nu_tst6 = np.array([])
y_nu_tst6 = np.array([])
xp_nu_tst6 = np.array([])
yp_nu_tst6 = np.array([])

E_nu_rf_tst6 = np.array([])
xp_nu_rf_tst6 = np.array([])
yp_nu_rf_tst6 = np.array([])

eH = eventHistory.eventHistory()
eH.outFile("Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst6.root")
eH.rootStructure()

for eventNum in range(n_tst6):
    #print("Event: ",eventNum)
    #print()

    #this is just needed for testing in order to reset the event history at the beginning of each loop
    #(otherwise locations that have not been filled this time but before will retain the earlier values)
    #eH = eventHistory.eventHistory()
    px = xppi[eventNum]*ppi[eventNum]
    py = yppi[eventNum]*ppi[eventNum]
    pz = ppi[eventNum]

    eH = nuTrfLineCmplx.GeneratePion(eventHist=eH,x=xpi[eventNum],y=ypi[eventNum],px=px,py=py,pz=pz,t=tpi[eventNum],weight=weight,runNum=runNum,eventNum=eventNum)

    pi_target = eH.findParticle('target')
    pi_prod = eH.findParticle('productionStraight')
    pi_dcy = eH.findParticle('pionDecay')
    mu = eH.findParticle('muonProduction')
    nu = eH.findParticle('piFlashNu')

    eH.fill()

    E_pi_trgt_tst6 = np.append(E_pi_trgt_tst6, pi_target.p()[0])

    if pi_target.traceSpace().s() != 0.:
        raise Exception("ERROR: s for pion at target is not 0 m!")


    ## Checks if pion hasn't decayed in the transfer line
    if pi_prod.run() != -1:

        t_pi_prod_tst6 = np.append(t_pi_prod_tst6, pi_prod.t())

        if pi_dcy.run() != -1 or mu.run() != -1 or nu.run() != -1:
            raise Exception("ERROR: Pion at production straight entrance generated even though pion decayed!")

        if pi_prod.traceSpace().s() != nuTrfLineCmplx.TrfLineCmplxLen():
            raise Exception("ERROR: s for pion at production straight entrance is not 67 m!")


    ## Checks if pion has decayed in the transfer line
    elif pi_dcy.run() != -1:

        n_dcy += 1

        z_pi_dcy_tst6 = np.append(z_pi_dcy_tst6, pi_dcy.traceSpace().z())
        t_pi_dcy_tst6 = np.append(t_pi_dcy_tst6, pi_dcy.t())

        E_mu_tst6 = np.append(E_mu_tst6, mu.p()[0])
        x_mu_tst6 = np.append(x_mu_tst6, mu.traceSpace().x())
        y_mu_tst6 = np.append(y_mu_tst6, mu.traceSpace().y())
        xp_mu_tst6 = np.append(xp_mu_tst6, mu.traceSpace().xp())
        yp_mu_tst6 = np.append(yp_mu_tst6, mu.traceSpace().yp())

        E_nu_tst6 = np.append(E_nu_tst6, nu.p()[0])
        x_nu_tst6 = np.append(x_nu_tst6, nu.traceSpace().x())
        y_nu_tst6 = np.append(y_nu_tst6, nu.traceSpace().y())
        xp_nu_tst6 = np.append(xp_nu_tst6, nu.traceSpace().xp())
        yp_nu_tst6 = np.append(yp_nu_tst6, nu.traceSpace().yp())

        P_mu_rf, P_numu_rf = nuTrfLineCmplx.BoostBack2RestFrame(PPi=ppi[eventNum],P_mu=mu.p(),P_numu=nu.p())

        #print("P_mu_rf: ",P_mu_rf)
        #print("P_numu_rf: ",P_numu_rf)
        #print()

        E_mu_rf_tst6 = np.append(E_mu_rf_tst6, P_mu_rf[0])
        xp_mu_rf_tst6 = np.append(xp_mu_rf_tst6, P_mu_rf[1][0]/P_mu_rf[1][2])
        yp_mu_rf_tst6 = np.append(yp_mu_rf_tst6, P_mu_rf[1][1]/P_mu_rf[1][2])

        E_nu_rf_tst6 = np.append(E_nu_rf_tst6, P_numu_rf[0])
        xp_nu_rf_tst6 = np.append(xp_nu_rf_tst6, P_numu_rf[1][0]/P_numu_rf[1][2])
        yp_nu_rf_tst6 = np.append(yp_nu_rf_tst6, P_numu_rf[1][1]/P_numu_rf[1][2])

        cosTheta = np.append(cosTheta,P_numu_rf[1][2]/P_numu_rf[0])
        phi_bf = np.arctan(P_numu_rf[1][1]/P_numu_rf[1][0])
        if (P_numu_rf[1][0] >= 0 and P_numu_rf[1][1] >= 0):
            phi = np.append(phi,phi_bf)
        elif (P_numu_rf[1][0] <= 0 and P_numu_rf[1][1] >= 0):
            phi = np.append(phi,phi_bf)
        elif (P_numu_rf[1][0] <= 0 and P_numu_rf[1][1] <= 0):
            phi = np.append(phi,phi_bf + np.pi/2)
        elif (P_numu_rf[1][0] >= 0 and P_numu_rf[1][1] <= 0):
            phi = np.append(phi,phi_bf - np.pi/2)

        if pi_prod.run() != -1:
            raise Exception("ERROR: Pion at production straight entrance generated even though pion decayed!")

        if mu.run() == -1 or nu.run() == -1:
            raise Exception("ERROR: Pion decayed but no daughter particle(s) generated.")

        if (pi_dcy.traceSpace().s() <= 0. or pi_dcy.traceSpace().s() >= nuTrfLineCmplx.TrfLineCmplxLen()):
            print("s = ",pi_dcy.traceSpace().s())
            raise Exception("ERROR: s for pion at decay point does not lie within range (0m,67m)!")

        #Check for energy conservation
        if mu.p()[0] + nu.p()[0] - pi_dcy.p()[0]  > 0.0001*pi_dcy.p()[0]:
            print(mu.p()[0] + nu.p()[0]," = ",pi_dcy.p()[0],"? Epi: ",Epi[eventNum])
            raise Exception("ERROR: Total energy is not conserved! (Maximum deviation allowed without raising exception: 0.01%)")

        #Check for 3-momentum conservation
        #print("Px: ",mu.p()[1][0]+nu.p()[1][0]," = ",pi_dcy.p()[1][0],"?")
        #print("Py: ",mu.p()[1][1]+nu.p()[1][1]," = ",pi_dcy.p()[1][1],"?")
        #if mu.p()[1][0] + nu.p()[1][0] - pi_dcy.p()[1][0] > 0.0001*pi_dcy.p()[1][0]:
        #    print(mu.p()[1][0]," + ",nu.p()[1][0]," = ",pi_dcy.p()[1][0],"?")
        #    raise Exception("ERROR: px is not conserved!")
        #
        #if mu.p()[1][1] + nu.p()[1][1] != pi_dcy.p()[1][1]:
        #    print(mu.p()[1][1]," + ",nu.p()[1][1]," = ",pi_dcy.p()[1][1],"?")
        #    raise Exception("ERROR: py is not conserved!")
        #
        #if mu.p()[1][2] + nu.p()[1][2] - pi_dcy.p()[1][2] > 0.0001*pi_dcy.p()[1][2]:
        #    print(mu.p()[1][2]+nu.p()[1][2]," = ",pi_dcy.p()[1][2],"?")
        #    raise Exception("ERROR: pz is not conserved! (Maximum deviation allowed without raising exception: 0.01%)")


print("    --> ",n_dcy," pions have decayed in the transfer line out of ")
print("        ",n_tst6," total pions at the target.")

eH.write()
eH.outFileClose()


##! Check boost back method
nuSTORMTrfLineCmplxTest = 7
print()
print("nuSTORMTrfLineCmplxTest:", nuSTORMTrfLineCmplxTest, " boost back method.")

PPi = 5.
p_mu = [0.01,0.02,0.1]
E_mu_tst7 = np.sqrt(p_mu[0]**2+p_mu[1]**2+p_mu[2]**2+(muCnst.mass()/1000.)**2)
P_mu = [E_mu_tst7,  p_mu]
p_numu = [-0.01,-0.02,-.1]
E_numu = np.sqrt(p_numu[0]**2+p_numu[1]**2+p_numu[2]**2)
P_numu = [E_numu,  p_numu]

P_mu_nS, P_numu_nS = nuTrfLineCmplx.Boost2nuSTORM(PPi,P_mu,P_numu)
P_mu_rf, P_numu_rf = nuTrfLineCmplx.BoostBack2RestFrame(PPi,P_mu_nS,P_numu_nS)

print("P_mu input: ",P_mu)
print("P_mu output after boosting back: ",P_mu_rf)
print("P_mu boosted to nuSTORM frame: ",P_mu_nS)
print()
print("P_numu input: ",P_numu)
print("P_numu output after boosting back: ",P_numu_rf)
print("P_numu boosted to nuSTORM frame: ",P_numu_nS)

'''
#! PLOTTING
#For nuSTORMTrfLineCmplxTest = 4
n, bins, patches = plt.hist(ppi, bins=50, color='y', range=(4.0,6.0))
plt.xlabel('Momentum (GeV)')
plt.ylabel('Entries')
plt.title('Momentum distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_p.pdf')
plt.close()

n, bins, patches = plt.hist(Epi, bins=50, color='y', range=(4.0,6.0))
plt.xlabel('Energy (GeV)')
plt.ylabel('Entries')
plt.title('Energy distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_E.pdf')
plt.close()

n, bins, patches = plt.hist(tpi, bins=int(t_bunch+t_spacing)*5, color='y', range=(0.,(t_bunch+t_spacing)*5*10**(-9)))
plt.xlabel('t (s)')
plt.ylabel('Entries')
plt.title('t distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_t.pdf')
plt.close()

n, bins, patches = plt.hist(-zpi-nuTrfLineCmplx.TrfLineCmplxLen(), bins=50, color='y', range=(-nuTrfLineCmplx.TrfLineCmplxLen()-max(zpi),-nuTrfLineCmplx.TrfLineCmplxLen()))
plt.xlabel('z (m)')
plt.ylabel('Entries')
plt.title('z distribution at t=0')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_z.pdf')
plt.close()

n, bins, patches = plt.hist(xpi, bins=50, color='y', range=(-0.16,0.16))
plt.xlabel('x (m)')
plt.ylabel('Entries')
plt.title('x distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_x.pdf')
plt.close()

n, bins, patches = plt.hist(ypi, bins=50, color='y', range=(-0.16,0.16))
plt.xlabel('y (m)')
plt.ylabel('Entries')
plt.title('y distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_y.pdf')
plt.close()

n, bins, patches = plt.hist(xppi, bins=50, color='y', range=(-0.0075,0.0075))
plt.xlabel('x^prime')
plt.ylabel('Entries')
plt.title('x^prime distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_xp.pdf')
plt.close()

n, bins, patches = plt.hist(yppi, bins=50, color='y', range=(-0.0075,0.0075))
plt.xlabel('y^prime')
plt.ylabel('Entries')
plt.title('y^prime distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_yp.pdf')
plt.close()

pltX = plt.hist2d(xpi,xppi, bins=200,range=[[-0.2,0.2],[-.008,.008]],norm=colors.LogNorm())
plt.colorbar()
plt.xlabel('x [m]')
plt.ylabel('x^prime')
plt.title('phase space distribution in x')
#plt.grid(alpha=0.)
#plt.show()
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_xPS_new.pdf')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_xPS_new.png')
plt.close()

n, bins, patches = plt.hist(t_dcy_pi, bins=50,color='y',range=(0.,0.7*10**(-5)),log=True) #,range=(0.,0.6*10**(-5))
plt.xlabel('$t_{decay}$ (s)')
plt.ylabel('Entries')
plt.title('Decay time distribution')
centers = bins[:-1] + np.diff(bins)/2
gamma_init = 100.
n0_init = n[0]
#l = 1./gamma
def special_exp(x,n0,gamma):
    y = n0*np.exp(-centers/(gamma*piCnst.lifetime()))
    return y
#plt.plot(centers, y, '-', color='b')
least_squares = LeastSquares(centers,n,yerror=np.sqrt(n),model=special_exp)
m = Minuit(least_squares,n0=n0_init,gamma=gamma_init)
m.migrad()
print(m.hesse())
plt.plot(centers,special_exp(centers,m.values[0],m.values[1]),color='b',label='fit')
#plt.show()
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_tdec.pdf')
plt.close()

n, bins, patches = plt.hist(s_dcy_pi, bins=50, color='y')
plt.xlabel('s decay (m)')
plt.ylabel('Entries')
plt.title('Decay point distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_sdec.pdf')
plt.close()

n, bins, patches = plt.hist(t_prod_pi-min(t_prod_pi), bins=50, color='y')
plt.xlabel('t (s) - '+str(min(t_prod_pi))+' s')
plt.ylabel('Entries')
plt.title('t distribution at production straight entrance')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/nuSTORMTrfLineCmplxTst_4_tprod.pdf')
plt.close()


'''
#For nuSTORMTrfLineCmplxTest = 5
dir = os.path.join(nuSIMPATH,r'Scratch/nuSTORMTrfLineCmplx/0_target')
dirExist = os.path.isdir(dir)
if dirExist == False:
    os.mkdir(dir)

n, bins, patches = plt.hist(E_pi_t, bins=50, color='y', range=(4.0,6.0))
plt.xlabel('Energy (GeV)')
plt.ylabel('Entries')
plt.title('Pion energy distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/0_target/nuSTORMTrfLineCmplxTst_5_piE.pdf')
plt.close()

n, bins, patches = plt.hist(t_pi_t, bins=int(t_bunch+t_spacing)*5, color='y', range=(0.,(t_bunch+t_spacing)*5*10**(-9)))
plt.xlabel('t (s)')
plt.ylabel('Entries')
plt.title('t distribution of pion at target')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/0_target/nuSTORMTrfLineCmplxTst_5_piT.pdf')
plt.close()

n, bins, patches = plt.hist(x_pi_t, bins=50, color='y', range=(-0.16,0.16))
plt.xlabel('x (m)')
plt.ylabel('Entries')
plt.title('x distribution of pion at target')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/0_target/nuSTORMTrfLineCmplxTst_5_piX.pdf')
plt.close()

n, bins, patches = plt.hist(y_pi_t, bins=50, color='y', range=(-0.16,0.16))
plt.xlabel('y (m)')
plt.ylabel('Entries')
plt.title('y distribution of pion at target')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/0_target/nuSTORMTrfLineCmplxTst_5_piY.pdf')
plt.close()

n, bins, patches = plt.hist(xp_pi_t, bins=50, color='y', range=(-0.0075,0.0075))
plt.xlabel('x^prime')
plt.ylabel('Entries')
plt.title('x^prime distribution of pion at target')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/0_target/nuSTORMTrfLineCmplxTst_5_piXp.pdf')
plt.close()

n, bins, patches = plt.hist(yp_pi_t, bins=50, color='y', range=(-0.0075,0.0075))
plt.xlabel('y^prime')
plt.ylabel('Entries')
plt.title('y^prime distribution of pion at target')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/0_target/nuSTORMTrfLineCmplxTst_5_piYp.pdf')
plt.close()


dir = os.path.join(nuSIMPATH,r'Scratch/nuSTORMTrfLineCmplx/1_productionStraight')
dirExist = os.path.isdir(dir)
if dirExist == False:
    os.mkdir(dir)

n, bins, patches = plt.hist(t_pi_prod-min(t_pi_prod), bins=50, color='y')
plt.xlabel('t (s) - '+str(min(t_pi_prod))+' s')
plt.ylabel('Entries')
plt.title('t distribution of pion at production straight entrance')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/1_productionStraight/nuSTORMTrfLineCmplxTst_5_piT.pdf')
plt.close()


dir = os.path.join(nuSIMPATH,r'Scratch/nuSTORMTrfLineCmplx/2_pionDecay')
dirExist = os.path.isdir(dir)
if dirExist == False:
    os.mkdir(dir)

n, bins, patches = plt.hist(t_pi_dcy, bins=50, color='y')
plt.xlabel('t (s)')
plt.ylabel('Entries')
plt.title('Pion decay time distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_5_piT.pdf')
plt.close()

n, bins, patches = plt.hist(s_pi_dcy, bins=50, color='y')
plt.xlabel('s decay (m)')
plt.ylabel('Entries')
plt.title('Pion decay point distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_5_piS.pdf')
plt.close()


n, bins, patches = plt.hist(E_mu, bins=50, color='y', range=(2.0,6.0))
plt.xlabel('Energy (GeV)')
plt.ylabel('Entries')
plt.title('Energy distribution of muon at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_5_muE.pdf')
plt.close()

n, bins, patches = plt.hist(x_mu, bins=50, color='y', range=(-0.16,0.16))
plt.xlabel('x (m)')
plt.ylabel('Entries')
plt.title('x distribution of muon at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_5_muX.pdf')
plt.close()

n, bins, patches = plt.hist(y_mu, bins=50, color='y', range=(-0.16,0.16))
plt.xlabel('y (m)')
plt.ylabel('Entries')
plt.title('y distribution of muon at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_5_muY.pdf')
plt.close()

n, bins, patches = plt.hist(xp_mu, bins=50, color='y', range=(-0.01,0.01))
plt.xlabel('x^prime')
plt.ylabel('Entries')
plt.title('x^prime distribution of muon at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_5_muXp.pdf')
plt.close()

n, bins, patches = plt.hist(yp_mu, bins=50, color='y', range=(-0.01,0.01))
plt.xlabel('y^prime')
plt.ylabel('Entries')
plt.title('y^prime distribution of muon at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_5_muYp.pdf')
plt.close()


n, bins, patches = plt.hist(E_nu, bins=50, color='y', range=(0.,2.5))
plt.xlabel('Energy (GeV)')
plt.ylabel('Entries')
plt.title('Energy distribution of PiFlashNu at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_5_nuE.pdf')
plt.close()

n, bins, patches = plt.hist(x_nu, bins=50, color='y', range=(-0.16,0.16))
plt.xlabel('x (m)')
plt.ylabel('Entries')
plt.title('x distribution of PiFlashNu at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_5_nuX.pdf')
plt.close()

n, bins, patches = plt.hist(y_nu, bins=50, color='y', range=(-0.16,0.16))
plt.xlabel('y (m)')
plt.ylabel('Entries')
plt.title('y distribution of PiFlashNu at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_5_nuY.pdf')
plt.close()

n, bins, patches = plt.hist(xp_nu, bins=50, color='y', range=(-0.1,0.1))
plt.xlabel('x^prime')
plt.ylabel('Entries')
plt.title('x^prime distribution of PiFlashNu at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_5_nuXp.pdf')
plt.close()

n, bins, patches = plt.hist(yp_nu, bins=50, color='y', range=(-0.1,0.1))
plt.xlabel('y^prime')
plt.ylabel('Entries')
plt.title('y^prime distribution of PiFlashNu at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_5_nuYp.pdf')
plt.close()


#For nuSTORMTrfLineCmplxTest = 6
n, bins, patches = plt.hist(E_pi_trgt_tst6, bins=50, color='y', range=(4.0,6.0))
plt.xlabel('Energy (GeV)')
plt.ylabel('Entries')
plt.title('Pion energy distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/0_target/nuSTORMTrfLineCmplxTst_6_piE.pdf')
plt.close()


n, bins, patches = plt.hist(t_pi_prod_tst6-min(t_pi_prod_tst6), bins=50, color='y')
plt.xlabel('t (s) - '+str(min(t_pi_prod_tst6))+' s')
plt.ylabel('Entries')
plt.title('t distribution of pion at production straight entrance')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/1_productionStraight/nuSTORMTrfLineCmplxTst_6_piT.pdf')
plt.close()


n, bins, patches = plt.hist(t_pi_dcy_tst6, bins=50, color='y')
plt.xlabel('t (s)')
plt.ylabel('Entries')
plt.title('Decay time distribution')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_piT.pdf')
plt.close()

n, bins, patches = plt.hist(z_pi_dcy_tst6, bins=50, color='y')
plt.xlabel('z decay (m)')
plt.ylabel('Entries')
plt.title('z distribution of pion decay')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_piZ.pdf')
plt.close()


n, bins, patches = plt.hist(E_mu_tst6, bins=50, color='y', range=(2.0,6.0))
plt.xlabel('Energy (GeV)')
plt.ylabel('Entries')
plt.title('Energy distribution of muon at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_muE.pdf')
plt.close()

n, bins, patches = plt.hist(x_mu_tst6, bins=50, color='y', range=(-0.16,0.16))
plt.xlabel('x (m)')
plt.ylabel('Entries')
plt.title('x distribution of muon at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_muX.pdf')
plt.close()

n, bins, patches = plt.hist(y_mu_tst6, bins=50, color='y', range=(-0.16,0.16))
plt.xlabel('y (m)')
plt.ylabel('Entries')
plt.title('y distribution of muon at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_muY.pdf')
plt.close()

n, bins, patches = plt.hist(xp_mu_tst6, bins=50, color='y', range=(-0.01,0.01))
plt.xlabel('x^prime')
plt.ylabel('Entries')
plt.title('x^prime distribution of muon at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_muXp.pdf')
plt.close()

n, bins, patches = plt.hist(yp_mu_tst6, bins=50, color='y', range=(-0.01,0.01))
plt.xlabel('y^prime')
plt.ylabel('Entries')
plt.title('y^prime distribution of muon at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_muYp.pdf')
plt.close()


n, bins, patches = plt.hist(E_nu_tst6, bins=50, color='y', range=(0.,2.5))
plt.xlabel('Energy (GeV)')
plt.ylabel('Entries')
plt.title('Energy distribution of PiFlashNu at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_nuE.pdf')
plt.close()

n, bins, patches = plt.hist(x_nu_tst6, bins=50, color='y', range=(-0.16,0.16))
plt.xlabel('x (m)')
plt.ylabel('Entries')
plt.title('x distribution of PiFlashNu at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_nuX.pdf')
plt.close()

n, bins, patches = plt.hist(y_nu_tst6, bins=50, color='y', range=(-0.16,0.16))
plt.xlabel('y (m)')
plt.ylabel('Entries')
plt.title('y distribution of PiFlashNu at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_nuY.pdf')
plt.close()

n, bins, patches = plt.hist(xp_nu_tst6, bins=50, color='y', range=(-0.1,0.1))
plt.xlabel('x^prime')
plt.ylabel('Entries')
plt.title('x^prime distribution of PiFlashNu at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_nuXp.pdf')
plt.close()

n, bins, patches = plt.hist(yp_nu_tst6, bins=50, color='y', range=(-0.1,0.1))
plt.xlabel('y^prime')
plt.ylabel('Entries')
plt.title('y^prime distribution of PiFlashNu at production point')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_nuYp.pdf')
plt.close()


n, bins, patches = plt.hist(E_mu_rf_tst6, bins=50, color='y', range=(0.1,0.12))
plt.xlabel('Energy (GeV)')
plt.ylabel('Entries')
plt.title('Energy distribution of muon at production point in pion restframe')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_muE_rf.pdf')
plt.close()

n, bins, patches = plt.hist(xp_mu_rf_tst6, bins=50, color='y', range=(-5.,5.))
plt.xlabel('x^prime')
plt.ylabel('Entries')
plt.title('x^prime distribution of muon at production point in pion restframe')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_muXp_rf.pdf')
plt.close()

n, bins, patches = plt.hist(yp_mu_rf_tst6, bins=50, color='y', range=(-5.,5.))
plt.xlabel('y^prime')
plt.ylabel('Entries')
plt.title('y^prime distribution of muon at production point in pion restframe')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_muYp_rf.pdf')
plt.close()


n, bins, patches = plt.hist(E_nu_rf_tst6, bins=50, color='y', range=(0.02,0.04))
plt.xlabel('Energy (GeV)')
plt.ylabel('Entries')
plt.title('Energy distribution of PiFlashNu at production point in pion restframe')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_nuE_rf.pdf')
plt.close()

n, bins, patches = plt.hist(xp_nu_rf_tst6, bins=50, color='y', range=(-5.,5.))
plt.xlabel('x^prime')
plt.ylabel('Entries')
plt.title('x^prime distribution of PiFlashNu at production point in pion restframe')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_nuXp_rf.pdf')
plt.close()

n, bins, patches = plt.hist(yp_nu_rf_tst6, bins=50, color='y', range=(-5.,5.))
plt.xlabel('y^prime')
plt.ylabel('Entries')
plt.title('y^prime distribution of PiFlashNu at production point in pion restframe')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_nuYp_rf.pdf')
plt.close()


n, bins, patches = plt.hist(cosTheta, bins=50, color='y', range=(-1.,1.))
plt.xlabel('$cos\Theta$')
plt.ylabel('Entries')
plt.title('$cos\Theta$ distribution for PiFlashNu')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_cosTheta_rf.pdf')
plt.close()


n, bins, patches = plt.hist(phi, bins=50, color='y') #, range=(-np.pi,np.pi)
plt.xlabel('$\phi$')
plt.ylabel('Entries')
plt.title('$\phi$ distribution for PiFlashNu')
plt.savefig('Scratch/nuSTORMTrfLineCmplx/2_pionDecay/nuSTORMTrfLineCmplxTst_6_phi_rf_tst.pdf')
plt.close()
'''

##! Complete:
print()
print("========  nuSTORMTrfLineCmplx: tests complete  ========")
