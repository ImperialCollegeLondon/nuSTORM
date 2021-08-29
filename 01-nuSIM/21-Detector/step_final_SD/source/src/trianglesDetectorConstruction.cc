//trianglesDetectorConstruction.cc
//Implementation of the trianglesDetectorConstruction class
//The definition of the class will be presented in (source directory of this project)/include/trianglesDetectorConstruction.hh
//
//
//
#include "trianglesDetectorConstruction.hh"
#include "trianglesScintSD.hh"
#include "stepSD.hh"
#include "trianglesMagneticField.hh"
#include "trianglesConstants.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4GenericTrap.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SystemOfUnits.hh"

#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4GenericMessenger.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

G4ThreadLocal trianglesMagneticField* trianglesDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* trianglesDetectorConstruction::fFieldMgr = 0;

//Initialize the constructor of the class.
trianglesDetectorConstruction::trianglesDetectorConstruction() : G4VUserDetectorConstruction(), fMessenger(nullptr), logicIron(nullptr){
	DefineCommands();
}

//Initialize the destructor of the class.
trianglesDetectorConstruction::~trianglesDetectorConstruction() {
	delete fMessenger;
}

// Implement the Construct() function which defines the geometry of the detector.
G4VPhysicalVolume* trianglesDetectorConstruction::Construct(){

	//Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();

	//Make sure the volumes do not overlap
	G4bool checkOverlaps = true;

	// Now we will define the geometry of the detector.
	// Two scintillators will form a Rhombus, and the Rhombuses will form a DetectorPlane (the words that start with capital letters are the names of the G4Solid).
	// Two DetectorPlanes with different orientations will be placed in a Module, and the Modules will form the WholeDetector.
	// The WholeDetector will be constructed in an orientation in which the beam direction would correspond to the y axis, because it is easier to construct the daughter volumes in that direction.
	// Finally the WholeDetector would be rotated and moved upon placement so that the beam direction would correspond to the z axis in the World.
	//
	// First, the size of the geometry of the scintillators.
	G4double triangle_base = 0.5*33. * mm;
	G4double base_center = 0.5 * triangle_base;
	G4double triangle_height = 0.5*17 * mm;
	G4double length = 4. * m;
	G4double height_hole = 0.5*8.5 * mm;

	// Next, the size of one detector plane.
	G4int num_repetition = kNofScintillators; // The number of rhombuses to create a plane
	G4double plane_sizeX = num_repetition * triangle_base; // The size of the detector plane.

	// And the size of the iron plates.
	G4double module_side = std::max(length, plane_sizeX+triangle_base);
	G4double iron_thickness = 1 * cm;
	G4Material* iron = nist->FindOrBuildMaterial("G4_Fe");

	// Then the size of one detector module
	G4double scintScint_distance = 0.0;
	G4double scintIron_distance = 0.0;
	G4double module_length = 4*scintScint_distance + 4*scintIron_distance + 6*triangle_height + 2*iron_thickness;
	// scintScint_distance and scintIron_distance are the thickness of the air between the detector planes.
	//
	//The size of the entire detector.
	G4int num_modules = kNofModules;
	G4double detector_length =  num_modules * module_length;

	// The number of scintillator strips facing one direction.
	fNoOfStrips = num_modules * 6 * num_repetition;
	//Finally define the size of the World.
	G4double world_sizeX = 1.5 * module_side;
	G4double world_sizeY = 1.5 * module_side;
	G4double world_sizeZ = 1.5 * detector_length;
	G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

	G4Box* solidWorld = new G4Box("World", 0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);
	// We will define all the solids first, and then put them all together at the end.
	//
	// First the entire detector.
	G4Box* solidWhole = new G4Box("Detector", 0.5*module_side, 0.5*detector_length, 0.5*module_side);

	//Next one detector module.
	G4Box* solidModule = new G4Box("Module", 0.5*module_side, 0.5*module_length, 0.5*module_side);

	// The iron plate.
	G4double diam_ironHole = 20*cm;
	G4Box* iron_box = new G4Box("Iron_Box", 0.5*module_side, 0.5*iron_thickness, 0.5*module_side);
	G4RotationMatrix* ironHole_rotate = new G4RotationMatrix;
	ironHole_rotate->rotateX(-(CLHEP::pi/2)*rad);
	G4Tubs* iron_hole = new G4Tubs("Iron_hole", 0.0, 0.5*diam_ironHole, 0.5*length, 0, CLHEP::twopi);
	G4SubtractionSolid* ironPlate = new G4SubtractionSolid("IronPlate", iron_box, iron_hole, ironHole_rotate, G4ThreeVector());

	// Then, we generate a detector plane using G4GenericTrap.
	std::vector<G4TwoVector> plane_verticies{G4TwoVector(-0.5*plane_sizeX,- 0.5*triangle_height), G4TwoVector(-0.5*plane_sizeX+base_center, 0.5*triangle_height), G4TwoVector(0.5*plane_sizeX+base_center, 0.5*triangle_height), G4TwoVector(0.5*plane_sizeX, -0.5*triangle_height), G4TwoVector(-0.5*plane_sizeX, -0.5*triangle_height), G4TwoVector(-0.5*plane_sizeX+base_center, 0.5*triangle_height), G4TwoVector(0.5*plane_sizeX+base_center, 0.5*triangle_height), G4TwoVector(0.5*plane_sizeX, -0.5*triangle_height)};
	G4GenericTrap* plane = new G4GenericTrap("Detector Plane", 0.5*length, plane_verticies);

	// Next, we generate a rhombus which includes two scintillators.
	// The verticies
	std::vector<G4TwoVector> rhombus_verticies{G4TwoVector(-base_center, -0.5*triangle_height), G4TwoVector(0.0, 0.5*triangle_height), G4TwoVector(triangle_base, 0.5*triangle_height), G4TwoVector(base_center, -0.5*triangle_height), G4TwoVector(-base_center, -0.5*triangle_height), G4TwoVector(0.0, 0.5*triangle_height), G4TwoVector(triangle_base, 0.5*triangle_height), G4TwoVector(base_center, -0.5*triangle_height)};
	// The solid
	G4GenericTrap* rhombus = new G4GenericTrap("Rhombus", 0.5*length, rhombus_verticies);

	// Then we generate a trianglular scintillator using G4TessellatedSolid.
	G4Material* scint_mat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

	// The verticies
	std::vector<G4TwoVector> scint_verticies{G4TwoVector(-base_center, -0.5*triangle_height), G4TwoVector(-base_center, -0.5*triangle_height), G4TwoVector(0.0, 0.5*triangle_height), G4TwoVector(base_center, -0.5*triangle_height), G4TwoVector(-base_center, -0.5*triangle_height), G4TwoVector(-base_center, -0.5*triangle_height), G4TwoVector(0.0, 0.5*triangle_height), G4TwoVector(base_center, -0.5*triangle_height)};
	//Then declare the solid.
	G4GenericTrap* scint_triangle = new G4GenericTrap("Scintillator_triangle", 0.5*length, scint_verticies);
	
	// Next we define the scintillator with the other rotation.
	// The verticies
	std::vector<G4TwoVector> scint2_verticies{G4TwoVector(base_center, -0.5*triangle_height), G4TwoVector(base_center, -0.5*triangle_height), G4TwoVector(0.0, 0.5*triangle_height), G4TwoVector(triangle_base, 0.5*triangle_height), G4TwoVector(base_center, -0.5*triangle_height), G4TwoVector(base_center, -0.5*triangle_height), G4TwoVector(0.0, 0.5*triangle_height), G4TwoVector(triangle_base, 0.5*triangle_height)};
	//Then declare the solid.
	G4GenericTrap* scint2_triangle = new G4GenericTrap("Scintillator2_triangle", 0.5*length, scint2_verticies);

	// Then we generate the cyllinders for the holes.
	G4double diam_scintHole = 2.6 * mm;
	G4Tubs* scintHole1 = new G4Tubs("Hole1", 0.0, 0.5*diam_scintHole, 0.5*length, 0, CLHEP::twopi);
	G4Tubs* scintHole2 = new G4Tubs("Hole2", 0.0, 0.5*diam_scintHole, 0.5*length, 0, CLHEP::twopi);

	// Transfer the holes...
	G4ThreeVector hole1_pos(0.0, -0.5*triangle_height+height_hole, 0.0);
	G4ThreeVector hole2_pos(base_center, 0.5*triangle_height-height_hole, 0.0);

	// Subtract the holes from the scintillators
	G4SubtractionSolid* scintillator1 = new G4SubtractionSolid("Scintillator1", scint_triangle, scintHole1, 0, hole1_pos);
	G4SubtractionSolid* scintillator2 = new G4SubtractionSolid("Scintillator2", scint2_triangle, scintHole2, 0, hole2_pos);

	// Generate the logical volumes
	logicWhole = new G4LogicalVolume(solidWhole, world_mat, "WholeDetector");
	logicModule = new G4LogicalVolume(solidModule, world_mat, "Module");
	logicPlane = new G4LogicalVolume(plane, world_mat, "DetectorPlane");
	logicIron = new G4LogicalVolume(ironPlate, iron, "IronPlate");
	logicRhombus = new G4LogicalVolume(rhombus, world_mat, "Rhombus");
	logicScint = new G4LogicalVolume(scintillator1, scint_mat, "Scintillator");
	logicScint2 = new G4LogicalVolume(scintillator2, scint_mat, "Scintillator2");
	logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");

	// Set the step limit in iron plate with magnetic field.
	G4UserLimits* userLimits = new G4UserLimits(1*m);
	logicIron->SetUserLimits(userLimits);

	// Place the scintillators in the Rhombus
	new G4PVPlacement(0, G4ThreeVector(), logicScint, "Scintillator", logicRhombus, false, 0, checkOverlaps);
	new G4PVPlacement(0, G4ThreeVector(), logicScint2, "Scintillator2", logicRhombus, false, 1, checkOverlaps);

	// Place the rhombuses in the detector plane.
	new G4PVReplica("Rhombuses", logicRhombus, logicPlane, kXAxis, num_repetition,triangle_base, 0);

	// Place the detector planes in the Module
	G4ThreeVector firstX(0.0, -0.5*module_length+0.5*triangle_height, 0.0);
	G4ThreeVector firstY(0.0, -0.5*module_length+1.5*triangle_height+scintScint_distance, 0.0);
	G4ThreeVector secondX(0.0, -0.5*module_length+2.5*triangle_height+2*scintScint_distance, 0.0);
	G4ThreeVector firstIron(0.0, -scintIron_distance-0.5*iron_thickness, 0.0);
	G4ThreeVector secondY(0.0, 0.5*triangle_height, 0.0);
	G4ThreeVector thirdX(0.0, 1.5*triangle_height+scintScint_distance, 0.0);
	G4ThreeVector thirdY(0.0, 2.5*triangle_height+2*scintScint_distance, 0.0);
	G4ThreeVector secondIron(0.0, 0.5*module_length-scintIron_distance-0.5*iron_thickness, 0.0);
	G4RotationMatrix* axistransfer = new G4RotationMatrix;
	axistransfer->rotateY((CLHEP::pi/2)*rad); // Apply this rotation matrix to the X scintillator planes.

	new G4PVPlacement(axistransfer, firstX, logicPlane, "DetectorPlaneX1", logicModule, false, 0, checkOverlaps);
	new G4PVPlacement(0, firstY, logicPlane, "DetectorPlaneY1", logicModule, false, 1, checkOverlaps);
	new G4PVPlacement(axistransfer, secondX, logicPlane, "DetectorPlaneX2", logicModule, false, 2, checkOverlaps);
	new G4PVPlacement(0, firstIron, logicIron, "Iron1", logicModule, false, 0, checkOverlaps);
	new G4PVPlacement(0, secondY, logicPlane, "DetectorPlaneY2", logicModule, false, 3, checkOverlaps);
	new G4PVPlacement(axistransfer, thirdX, logicPlane, "DetectorPlaneX3", logicModule, false, 4, checkOverlaps);
	new G4PVPlacement(0, thirdY, logicPlane, "DetectorPlaneY3", logicModule, false, 5, checkOverlaps);
	new G4PVPlacement(0, secondIron, logicIron, "Iron2", logicModule, false, 1, checkOverlaps);
	
	// Place the Modules in the Detector.
	new G4PVReplica("Modules", logicModule, logicWhole, kYAxis, num_modules, module_length, 0);

	// Place the Detector in the world.
		// Rotate the Detector
		G4RotationMatrix* whole_world = new G4RotationMatrix;
		whole_world->rotateX(-(CLHEP::pi/2)*rad);
	new G4PVPlacement(whole_world, G4ThreeVector(), logicWhole, "Whole", logicWorld, false, 0, checkOverlaps);
	
	// Finally, place the World in the appropriate position.
	G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);

	// Set the visualization attributes.
	G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	visAttributes->SetVisibility(false);
	logicWorld->SetVisAttributes(visAttributes);
	fVisAttributes.push_back(visAttributes); // Add the visualization attributes of the World to the end of the array fVisAttributes.

	visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	visAttributes->SetVisibility(false);
	logicWhole->SetVisAttributes(visAttributes);
	fVisAttributes.push_back(visAttributes);

	visAttributes = new G4VisAttributes(G4Colour(0.9, 0.9, 0.9));
	logicIron->SetVisAttributes(visAttributes);
	fVisAttributes.push_back(visAttributes);

	visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 0.3));
	logicPlane->SetVisAttributes(visAttributes);
	fVisAttributes.push_back(visAttributes);

	visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	visAttributes->SetVisibility(false);
	logicModule->SetVisAttributes(visAttributes);
	fVisAttributes.push_back(visAttributes);

	visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	visAttributes->SetVisibility(false);
	logicRhombus->SetVisAttributes(visAttributes);
	logicScint->SetVisAttributes(visAttributes);
	logicScint2->SetVisAttributes(visAttributes);
	fVisAttributes.push_back(visAttributes);


	return physWorld;

}

void trianglesDetectorConstruction::ConstructSDandField(){
	// The Sensitive detectors.
	// The scintillator strips are specified by the copy numbers of the module (0~num_modules-1), plane (0~5), rhombus (0~num_repetition-1) and then the scintillator itself (0,1).
	auto sdManager = G4SDManager::GetSDMpointer();
	G4String SDname;
	auto scintillators = new trianglesScintSD(SDname="/scintillators", "scintColl", fNoOfStrips);
	sdManager->AddNewDetector(scintillators);
	SetSensitiveDetector("Scintillator", scintillators);

	auto scintillators2 = new trianglesScintSD("/scintillators2", "scint2Coll", fNoOfStrips);
	sdManager->AddNewDetector(scintillators2);
	SetSensitiveDetector("Scintillator2", scintillators2);

	auto scint_steps = new stepSD("/scint_steps", "scint_stepsColl");
	sdManager->AddNewDetector(scint_steps);
	SetSensitiveDetector("Scintillator", scint_steps);

	auto scint2_steps = new stepSD("/scint2_steps", "scint2_stepsColl");
	sdManager->AddNewDetector(scint2_steps);
	SetSensitiveDetector("Scintillator2", scint2_steps);

	// The magnetic field. (Sensitive detectors might be added in future updates)
	fMagneticField = new trianglesMagneticField();
	fFieldMgr = new G4FieldManager();
	fFieldMgr->SetDetectorField(fMagneticField);
	fFieldMgr->CreateChordFinder(fMagneticField);
	G4bool forceToAllDaughters = true;
	logicIron->SetFieldManager(fFieldMgr, forceToAllDaughters);

}

void trianglesDetectorConstruction::DefineCommands(){
	fMessenger = new G4GenericMessenger(this, "/triangles/detector/", "Detector control");

}
