// trianglesDetectorConstruction.hh
// Brief defenition of the trianglesDetectorConstruction class.
// The inplementation should be found in (source directory of this project)/src/trianglesDetectorConstruction.cc
//
#ifndef trianglesDetectorConstruction_h
#define trianglesDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4FieldManager.hh"

class trianglesMagneticField;

class G4VPhysicalVolume;
//class G4LogicalVolume;
class G4VisAttributes;
class G4GenericMessenger;

// Detector construction
class trianglesDetectorConstruction : public G4VUserDetectorConstruction{

	public:
		trianglesDetectorConstruction();
		virtual ~trianglesDetectorConstruction();

		virtual G4VPhysicalVolume* Construct();

		virtual void ConstructSDandField();

	private:
		void DefineCommands();

		G4GenericMessenger* fMessenger;

		static G4ThreadLocal trianglesMagneticField* fMagneticField;
		static G4ThreadLocal G4FieldManager* fFieldMgr;

		G4LogicalVolume* logicWhole;
		G4LogicalVolume* logicModule;
		G4LogicalVolume* logicPlane;
		G4LogicalVolume* logicIron;
		G4LogicalVolume* logicRhombus;
		G4LogicalVolume* logicScint;
		G4LogicalVolume* logicScint2;
		G4LogicalVolume* logicWorld;

		G4int fNoOfStrips; // Number of scintillator strips.
		std::vector<G4VisAttributes*> fVisAttributes;
};

#endif
