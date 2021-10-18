#!/usr/bin/env python3

"""
Python script to create the target and horn geometry
for BDSIM and Fluka using pyg4ometry.

Version history:
----------------------------------------------
 1.0: 13Oct21: First version using Inconel target

"""

# nuSTORM target and horn geometry
import pyg4ometry
import pyg4ometry.geant4 as g4
import pyg4ometry.gdml as gdml
import pyg4ometry.fluka as fluka
import pyg4ometry.convert as convert

import math
import os

# Math constants
pi = math.acos(-1.0)
twoPi = 2.0*pi
# Zero rotation
zeroRot = [0.0, 0.0, 0.0]
# Geant4 uses mm, Fluka cm
mm2cm = 0.1

# Use mm dimensions throughout: default Geant4 position unit.
# Otherwise, add 'unit' argument at end of volumes & position definitions
# e.g. Tubs(..., gReg, 'cm') and [0.0, 0.0, 100.0, 'cm']

# Make the visualisation of cylinders a little smoother
pyg4ometry.config.SolidDefaults.Tubs.nslice = 32

class hornPars():

    def __init__(self, pol = 1):

        # Current (kA)
        self.I = 230.0
        # Current polarity: +1 or -1
        self.pol = pol
        # Fluka current parameter file: fixed in Fluka/src/user/magfld.f
        self.flukaBParFile = 'FieldPars.dat'
        
        # Radii
        self.rOut = 100.0
        self.rUps = 85.0
        self.rIn  = 20.0 # 10.0
        self.rDwn = 50.0

        # Lengths
        self.LTot = 2200.0
        self.LUps = 250.0
        self.LMid = 850.0
        self.LDwn = self.LTot - self.LUps - self.LMid

        # Upstream polycone points
        nUps = 20
        dPhi = 0.5*pi/(nUps*1.0)
        dzUpsIC = self.LUps/(nUps*1.0)
        drUpsIC = (self.rUps - self.rIn)
        self.zUpsIC = []
        self.rInUpsIC = []
        self.rOutUpsIC = []
        for i in range(nUps+1):
            z = dzUpsIC*i
            r = self.rIn + drUpsIC*math.cos(dPhi*i)
            self.zUpsIC.append(z)
            self.rInUpsIC.append(0.0)
            self.rOutUpsIC.append(r)
            
        #print('zUpsIC = {0}'.format(self.zUpsIC))
        #print('rUpsIC = {0}'.format(self.rOutUpsIC))

        # Downstream polycone points
        nDwn = 20
        dPhi = 0.5*pi/(nDwn*1.0)
        dzDwnIC = self.LDwn/(nDwn*1.0)
        drDwnIC = (self.rDwn - self.rIn)
        self.zDwnIC = []
        self.rInDwnIC = []
        self.rOutDwnIC = []
        for i in range(nDwn+1):
            z = dzDwnIC*i + self.LUps + self.LMid
            r = self.rIn + drDwnIC*math.sin(dPhi*i)
            self.zDwnIC.append(z)
            self.rInDwnIC.append(0.0)
            self.rOutDwnIC.append(r)
            
        #print('zDwnIC = {0}'.format(self.zDwnIC))
        #print('rDwnIC = {0}'.format(self.rOutDwnIC))

        # Conductor thickness: 2.5 mm
        self.thick = 2.5

        # Polycone points for the upstream inner Ar gas field region
        self.zUpsAr = []
        self.rInUpsAr = []
        self.rOutUpsAr = []
        dzUpsAr = (self.LUps - self.thick)/(nUps*1.0)
        drUpsAr = (self.rUps - self.rIn)
        for i in range(nUps+1):
            z = dzUpsAr*i
            r = self.rIn + drUpsIC*math.cos(dPhi*i) + self.thick
            self.zUpsAr.append(z)
            self.rInUpsAr.append(0.0)
            self.rOutUpsAr.append(r)

        #print('zUpsAr = {0}'.format(self.zUpsAr))
        #print('rInUpsAr = {0}'.format(self.rInUpsAr))
        #print('rOutUpsAr = {0}'.format(self.rOutUpsAr))

        # Polycone points for the downstream inner Ar gas field region
        self.zDwnAr = []
        self.rInDwnAr = []
        self.rOutDwnAr = []
        dzDwnAr = (self.LDwn - self.thick)/(nDwn*1.0)
        drDwnAr = (self.rDwn - self.rIn)
        for i in range(nDwn+1):
            z = dzDwnAr*i + self.LUps + self.LMid
            r = self.rIn + drDwnAr*math.sin(dPhi*i) + self.thick
            self.zDwnAr.append(z)
            self.rInDwnAr.append(0.0)
            self.rOutDwnAr.append(r)

        #print('zDwnAr = {0}'.format(self.zDwnAr))
        #print('rInDwnAr = {0}'.format(self.rInDwnAr))
        #print('rOutDwnAr = {0}'.format(self.rOutDwnAr))


class targetPars():

    def __init__(self, hornPars):

        # Cylinder: r = 6.3 mm, L = 46 cm
        self.r = 6.3
        self.L = 460.0
        # Start of target from z = 0 (mm)
        # nuSTORM thesis 550.0 mm
        self.zStart = hornPars.LUps

        self.material = 'Inconel'

        
class beamPars():

    def __init__(self, beamE):

        # Beam parameters: KE (GeV), sigma (mm), z0 (mm)
        self.E = beamE
        self.sigma = 2.67
        self.FWHM = 2.0*math.sqrt(2.0*math.log(2.0))*self.sigma
        self.z0 = -100.0

        # Number of events
        self.N = 1000

        
class allParameters():

    def __init__(self, hornPars, targetPars, beamPars):

        self.hornPars = hornPars
        self.targetPars = targetPars
        self.beamPars = beamPars

        
class rzUsrbin():

    def __init__(self, name, part, unit, rMin, rMax, dr, zMin, zMax, dz):

        self.name = name
        self.part = part
        # Negative unit for binary
        self.unit = unit
        # r bins
        self.rMin = rMin
        self.rMax = rMax
        self.dr = dr
        self.Nr = int((rMax - rMin)/(dr*1.0) + 0.5)
        # z bins
        self.zMin = zMin
        self.zMax = zMax
        self.dz = dz
        self.Nz = int((zMax - zMin)/(dz*1.0) + 0.5)
        # 1 phi bin
        self.x = 0.0
        self.y = 0.0
        self.Nphi = 1
        #print('Nr = {0}, Nz = {1}'.format(self.Nr, self.Nz))


def setMaterials(gReg):

    # Using NIST material compound definitions where possible.
    # Must use Fluka capital names, so that their neutron cross-sections are automatically used.
    # Otherwise, we would need to associate materials with appropriate LOW-MAT definitions.
    # Fluka conversion adds unnecessary lines, so we remove these in addFlukaRunInfo().
    # Fluka cross-sections: http://www.fluka.org/content/manuals/online/10.4.1.2.html
    # These materials should also work for Geant4
    
    air = g4.MaterialCompound('AIR', 1.205e-3, 2, gReg)
    N = g4.MaterialSingleElement('NITROGEN', 7, 14.0, 1.165e-3, gReg)
    O = g4.MaterialSingleElement('OXYGEN', 8, 16.0, 1.332e-3, gReg)
    air.add_material(N, 0.76)
    air.add_material(O, 0.24)

    argon = g4.MaterialSingleElement('ARGON', 18, 38.0, 1.662e-3, gReg)
    
    Al = g4.MaterialSingleElement('ALUMINIUM', 13, 27, 2.7, gReg)

    # Inconel-625: http://www.goodfellow.com/A/Inconel-625-Corrosion-Resistant-Alloy.html
    inconel = g4.MaterialCompound('INCONEL', 8.44, 6, gReg)
    Ni = g4.MaterialSingleElement('NICKEL', 28, 61.0, 8.90, gReg)
    Cr = g4.MaterialSingleElement('CHROMIUM', 24, 52.0, 7.18, gReg)
    Mo = g4.MaterialSingleElement('MOLYBDENUM', 42, 96.0, 10.22, gReg)
    Fe = g4.MaterialSingleElement('IRON', 26, 56.0, 7.874, gReg)
    Nb = g4.MaterialSingleElement('NIOBIUM', 41, 93.0, 8.57, gReg)
    Ti = g4.MaterialSingleElement('TITANIUM', 22, 48.0, 4.54, gReg)    
    inconel.add_material(Ni, 0.61)
    inconel.add_material(Cr, 0.22)
    inconel.add_material(Mo, 0.09)
    inconel.add_material(Fe, 0.05)
    # Add 3% of Nb & Ti to get 100% total
    inconel.add_material(Nb, 0.015)
    inconel.add_material(Ti, 0.015)


def run(beamE, pol):

    print('Running nuSTORMTarget.py: beamE = {0} GeV and horn polarity = {1}'.format(beamE, pol))

    # Geant4 registry to store gdml data
    gReg = g4.Registry()

    # Geometry parameters
    hP = hornPars(pol)
    tP = targetPars(hP)

    # Materials
    setMaterials(gReg)
        
    # World volume and logical
    worldLxy = 4.0*hP.rOut
    worldLz  = 4.0*hP.LTot
    world  = g4.solid.Box('World', worldLxy, worldLxy, worldLz, gReg)
    worldL = g4.LogicalVolume(world, 'AIR', 'WorldL', gReg)
    gReg.setWorld(worldL.name)

    # Create horn and target
    createHorn(hP, gReg)
    createTarget(tP, gReg)
    # Define aperture volume for tallying particles
    # leaving target & horn downstream area
    defineAperture(hP, gReg)
    
    # Check for overlaps
    #worldL.checkOverlaps(recursive=True)

    # Visualise geometry
    v = pyg4ometry.visualisation.VtkViewer()
    v.addLogicalVolume(worldL)
    v.addAxes(length=10)
    v.view()

    # Write gdml file
    w = gdml.Writer()
    w.addDetector(gReg)
    w.write('nuSTORMTarget.gdml')

    # Convert geometry to Fluka format. First write to a file
    # only containing the geometry, then add other run options.
    flukaFile = 'nuSTORMTarget.inp'

    # Geometry
    flukaGeomFile = flukaFile.replace('.inp', 'Geom.inp')
    fReg = convert.geant4Reg2FlukaReg(gReg)
    wFluka = fluka.Writer()
    wFluka.addDetector(fReg)
    print('Writing {0}'.format(flukaGeomFile))
    wFluka.write(flukaGeomFile)

    # Set B field for horn region
    flukaGeomBFile = flukaGeomFile.replace('.inp', '_B.inp')
    BParFile = setFlukaBFields(flukaGeomFile, flukaGeomBFile, hP)
    
    # Add p beam, B field & physics run options to Fluka file
    bP = beamPars(beamE)
    allPars = allParameters(hP, tP, bP)
    addFlukaRunInfo(flukaFile, flukaGeomBFile, allPars)
        
    # Create Flair template file
    extent = worldL.extent(includeBoundingSolid=True)
    flairFile = flukaFile.replace('.inp', '.flair')
    flair = fluka.Flair(flukaFile, extent)
    print('Writing {0}'.format(flairFile))
    flair.write(flairFile)
    
    # BDSIM job template
    w.writeGmadTester('nuSTORMTarget.gmad', 'nuSTORMTarget.gdml')


def createHorn(hP, gReg):
    
    print('Creating Horn')

    worldL = gReg.getWorldVolume()
    
    # Horn main cylinder volume, defining inner and outer conductor ranges
    hornConCyl = g4.solid.Tubs('hornConCyl', hP.rIn, hP.rOut, hP.LTot, 0.0, twoPi, gReg)

    # Upstream inner conductor polycone
    hornUpsIC = g4.solid.Polycone('hornUpsIC', 0.0, twoPi, hP.zUpsIC, hP.rInUpsIC, hP.rOutUpsIC, gReg)

    # Downstream inner conductor polycone
    hornDwnIC = g4.solid.Polycone('hornDwnIC', 0.0, twoPi, hP.zDwnIC, hP.rInDwnIC, hP.rOutDwnIC, gReg)

    # Subtract upstream polycone from main cylinder
    posUpsIC = [0.0, 0.0, -0.5*hP.LTot]
    hornConA  = g4.solid.Subtraction('hornConA', hornConCyl, hornUpsIC, [zeroRot, posUpsIC], gReg)

    # Subtract downstream polycone from previous shape
    posDwnIC = [0.0, 0.0, -0.5*hP.LTot]
    hornCon  = g4.solid.Subtraction('hornCon', hornConA, hornDwnIC, [zeroRot, posDwnIC], gReg)

    # Logical and physical volumes
    hornConL = g4.LogicalVolume(hornCon, 'ALUMINIUM', 'hornConL', gReg)
    hornConPos = [0.0, 0.0, 0.5*hP.LTot]
    hornConPhys = g4.PhysicalVolume(zeroRot, hornConPos, hornConL, 'hornConPhys', worldL, gReg)

    # Horn Ar gas main cylinder volume
    hornArCyl = g4.solid.Tubs('hornArCyl', hP.rIn + hP.thick, hP.rOut - hP.thick,
                              hP.LTot - 2.0*hP.thick, 0.0, twoPi, gReg)

    # Upstream inner gas polycone
    hornUpsAr = g4.solid.Polycone('hornUpsAr', 0.0, twoPi, hP.zUpsAr, hP.rInUpsAr, hP.rOutUpsAr, gReg)

    # Downstream inner gas polycone
    hornDwnAr = g4.solid.Polycone('hornDwnAr', 0.0, twoPi, hP.zDwnAr, hP.rInDwnAr, hP.rOutDwnAr, gReg)

    # Subtract upstream polycone from main cylinder
    posUpsAr = [0.0, 0.0, -0.5*hP.LTot + hP.thick]
    hornArA  = g4.solid.Subtraction('hornArA', hornArCyl, hornUpsAr, [zeroRot, posUpsAr], gReg)

    # Subtract downstream polycone from previous shape
    posDwnAr = [0.0, 0.0, -0.5*hP.LTot + hP.thick]
    hornAr  = g4.solid.Subtraction('hornAr', hornArA, hornDwnAr, [zeroRot, posDwnAr], gReg)

    # Logical and physical volumes
    hornArL = g4.LogicalVolume(hornAr, 'ARGON', 'hornArL', gReg)
    hornArPos = [0.0, 0.0, 0.0]
    hornConPhys = g4.PhysicalVolume(zeroRot, hornArPos, hornArL, 'hornArPhys', hornConL, gReg)

    
def createTarget(tP, gReg):

    print('Creating target')

    worldL = gReg.getWorldVolume()

    # Target cylinder
    targetCyl = g4.solid.Tubs('targetCyl', 0.0, tP.r, tP.L, 0.0, twoPi, gReg)

    # Logical and physical volumes
    targetCylL = g4.LogicalVolume(targetCyl, 'INCONEL', 'targetCylL', gReg)
    targetPos  = [0.0, 0.0, tP.zStart + 0.5*tP.L]
    targetPhys = g4.PhysicalVolume(zeroRot, targetPos, targetCylL, 'targetPhys', worldL, gReg)


def defineAperture(hP, gReg):

    print('Defining 10 mm long tracking aperture')

    worldL = gReg.getWorldVolume()

    aL = 2.0
    # Cylindrical aperture shape
    aperture = g4.solid.Tubs('AIR', 0.0, hP.rOut, aL, 0.0, twoPi, gReg)

    # Logical and physical volumes
    apertureL = g4.LogicalVolume(aperture, 'AIR', 'apertureL', gReg)
    aperturePos = [0.0, 0.0, hP.LTot + 0.5*aL]
    aperturePhys = g4.PhysicalVolume(zeroRot, aperturePos, apertureL, 'aperturePhys', worldL, gReg)
    
    
def setFlukaBFields(flukaGeomFile, flukaGeomBFile, hP):

    print('Setting Fluka B fields for {0}'.format(flukaGeomFile))
    inFile = open(flukaGeomFile, 'r')
    outFile = open(flukaGeomBFile, 'w')
    
    for line in inFile:

        # Only need B field inside horn argon gas region
        if 'ASSIGNMA' in line and 'ARGON' in line:
            BLine = line.rstrip()
            BLine = BLine.replace(' ', ', ')
            BLine = BLine + ', , , 1.\n'
            outFile.write(BLine)
        else:
            outFile.write(line)

    outFile.close()
    inFile.close()

    # Set horn current parameters
    parFile = open(hP.flukaBParFile, 'w')
    parFile.write('{0} {1}\n'.format(hP.I, hP.pol))
    parFile.close()

        
def addFlukaRunInfo(flukaFile, geomFile, allPars):

    print('Defining {0} with geometry from {1}'.format(flukaFile, geomFile))
    
    outFile = open(flukaFile, 'w')

    # Title and defaults
    outFile.write('TITLE\n')
    outFile.write('nuSTORM target simulation\n')
    outFile.write('DEFAULTS {0:>69}\n'.format('PRECISIO'))

    # Beam energy, size and particle type
    bP = allPars.beamPars
    line = '{0:<10}{1:>10.1f}{2:>10}{3:>10}{4:>10.3f}{5:>10.3f}{6:>16}\n'.format('BEAM', -bP.E, 0.0, 0.0,
                                                                                 -bP.FWHM*mm2cm,
                                                                                 -bP.FWHM*mm2cm,
                                                                                 'PROTON')
    outFile.write(line)

    # Beam posiiton and direction
    line = '{0:<10}{1:>10.1f}{2:>10}{3:>10}{4:>10.1f}{5:>10.1f}\n'.format('BEAMPOS', 0.0, 0.0,
                                                                          bP.z0*mm2cm, 0.0, 0.0)
    outFile.write(line)

    # Set magnetic field tracking parameters
    line = '{0:<10}{1:>10.1f}{2:>10.2f}{3:>10.2f}{4:>10.1f}{5:>10.1f}{6:>10.1f}\n'.format('MGNFIELD', 30.0, 0.05,
                                                                                          0.05, 0.0, 0.0, 0.0)
    outFile.write(line)

    # Write geometry
    NLines = sum(1 for line in open(geomFile, 'r'))

    matList = []
    
    with open(geomFile) as geom:
        for lineNo, line in enumerate(geom, 1):
            # Ignore last END statement
            if lineNo < NLines:
                # Remove extraneous material definitions
                # which get added by the Fluka conversion
                (ok, newLine) = checkMatLine(line, matList)
                if ok == True:
                    outFile.write(newLine)
        
    # Physics: PEANUT model energy limits
    EMax = 1000.0
    line = 'PHYSICS, {0:.1f}, {1:.1f}, {2:.1f}, {3:.1f}, {4:.1f}, {5:.1f}, ' \
           'PEATHRES\n'.format(EMax, EMax, EMax, EMax, EMax, EMax)
    outFile.write(line)

    # Enable evaporation
    outFile.write('PHYSICS, 3.0, , , , , , EVAPORAT\n')
    
    # Thresholds
    outFile.write('THRESHOL, , , 0.0, 0.0, , 0.0\n')

    # Enable coalescence, which requires linking with DPMJET
    outFile.write('PHYSICS, 1.0, , , , , , COALESCE\n')

    # Required for deuterons only
    outFile.write('PHYSICS, 1.0, 0.001, 0.15, 2.0, 2.0, 3.0, IONSPLIT\n')

    # Transport ions (from radio-nuclear interactions)
    outFile.write('IONTRANS, HEAVYION\n')

    # Score energy deposition in all regions for each event
    outFile.write('SCORE, ENERGY\n')

    # Dump regional event scoring info to output file
    outFile.write('EVENTDAT, 60.0, , , , , , evt.txt\n')
    
    # Usrbin histograms
    hP = allPars.hornPars
    tP = allPars.targetPars
    
    # Energy density: all, target
    EAll = rzUsrbin('EALL', 'ENERGY', -21.0, 0.0, hP.rOut*mm2cm+10.0, 0.1,
                    bP.z0*mm2cm, hP.LTot*mm2cm+10.0, 1.0)

    ETgt = rzUsrbin('ETGT', 'ENERGY', -22.0, 0.0, tP.r*mm2cm, 0.01,
                    tP.zStart*mm2cm, (tP.L+tP.zStart)*mm2cm, 1.0)

    # List of all usrbins
    ubinList = [EAll, ETgt]
    for iB,u in enumerate(ubinList):
    
        line = 'USRBIN, 11.0, {0}, {1:.1f}, {2:.2f}, {3:.2f}, ' \
               '{4:.2f}, {5}\n'.format(u.part, u.unit, u.rMax, u.y, u.zMax, u.name)        
        outFile.write(line)

        line = 'USRBIN, {0:.2f}, {1:.2f}, {2:.2f}, {3:.2f}, {4:.2f}, ' \
               '{5:.2f} &\n'.format(u.rMin, u.x, u.zMin, u.Nr, u.Nphi, u.Nz)
        outFile.write(line)

    # Dump particle info passing the tracking aperture via Fluka/src/user/mgdraw.f.
    # Userdump option = 3 only uses BXDRAW output, since EEDRAW does nothing and we
    # comment out the ENDRAW code. The file "particles" (unit 42) should be empty,
    # since we only write BXDRAW output that is directed to "output.txt" (unit 41)
    outFile.write('USERDUMP, 100.0, 42.0, 3.0, 1.0, , , particles\n')
    outFile.write('OPEN, 41.0, , , , , , UNKNOWN\n')
    outFile.write('output.txt\n')
    
    # Random seed
    outFile.write('RANDOMIZ, 1.0, 1.0\n')

    # Number of primaries
    line = 'START, {0}.\n'.format(bP.N)
    outFile.write(line)

    # Finalise
    outFile.write('STOP\n')
    outFile.close()
    

def checkMatLine(line, matList):

    isOK = True

    # Replace ALUMINIU with ALUMINUM
    newLine = line.replace('ALUMINIU', 'ALUMINUM')
    sLine = newLine.split()
    N = len(sLine)
    
    if N == 8:
        if sLine[0] == 'MATERIAL,':
            # Remove extraneous "0" numbers
            matName = sLine[-1]
            matName = matName.replace('0', '')
            if matName in matList:
                # Already have this material
                isOK = False
            else:
                matList.append(matName)
                sLine[-1] = matName
                # Update line
                newLine = ''.join(sLine) + '\n'

        elif sLine[0] == 'COMPOUND,':
            # Remove extraneous "0" numbers from names
            matName1 = sLine[2].replace('0', '')
            matName2 = sLine[4].replace('0', '')
            matName3 = sLine[6].replace('0', '')
            matName4 = sLine[7].replace('0', '')
            newLine = 'COMPOUND, {0} {1} {2} {3} {4} {5} {6} \n'.format(sLine[1], matName1, sLine[3],
                                                                        matName2, sLine[5], matName3,
                                                                        matName4)            
    return (isOK, newLine)


if __name__ == '__main__':

    # Proton KE (GeV) and horn field polarity (+1 or -1)
    beamE = 100.0
    pol = 1.0

    nArg = len(os.sys.argv)
    if nArg > 1:
        beamE = float(os.sys.argv[1])
    if nArg > 2:
        pol = float(os.sys.argv[2])
        
    run(beamE, pol)
