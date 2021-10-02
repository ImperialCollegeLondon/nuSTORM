// trianglesMagneticField.cc
// Implementation of the trianglesMagneticField class. See the corresponding header files for the definition, purpose, etc. of this class.
//
#include "trianglesMagneticField.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

trianglesMagneticField::trianglesMagneticField() : G4MagneticField(), fMessenger(nullptr), fB(1.5*tesla){
	// Define the commands for this class upon construction of the instance.
	DefineCommands();
}

trianglesMagneticField::~trianglesMagneticField(){
	delete fMessenger;
}

void trianglesMagneticField::GetFieldValue(const G4double [4], double *bField) const{
	// bfield is the value of the field at each spacetime point.
	bField[0] = 0.0;
	bField[1] = fB;
	bField[2] = 0.0;
}

void trianglesMagneticField::DefineCommands(){
	// Define /triangles/field command directory using the generic messenger class
	fMessenger = new G4GenericMessenger(this, "/triangles/field/", "Field control");

	// fieldValue command
	auto& valueCmd = fMessenger->DeclareMethodWithUnit("value", "tesla", &trianglesMagneticField::SetField, "Set field strength.");
	// (above) Declares a method with the name "value", which uses tesla as the default unit. 
	// When you type /triangles/field/set 4, the method will run trianglesMagneticField::SetField with val = 4. 
	// "Set field Strength." is the explanation of the command you can see when you execute ls /B5/field
	
	valueCmd.SetParameterName("field", true);
	valueCmd.SetDefaultValue("1.5");
}
