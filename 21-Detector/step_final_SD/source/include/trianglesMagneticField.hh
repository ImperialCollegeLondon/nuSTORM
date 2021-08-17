// trianglesMagneticField.hh
// Definition of the trianglesMagneticField class.
// This class defines the magnetic fields that will be invoked to the iron plates.
//
#ifndef B5MagneticField_H
#define B5MagneticField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"

class G4GenericMessenger;

class trianglesMagneticField : public G4MagneticField{
	public:
		trianglesMagneticField();
		virtual ~trianglesMagneticField();

		// Get the value of the magnetic field at a given spacetime point.
		virtual void GetFieldValue(const G4double point[4], double* bField) const;

		// Set the field strength to a value of your choice (the imput from a command?)
		void SetField(G4double val) { fB = val;}
		G4double GetField() const { return fB;}
	private:
		void DefineCommands();

		G4GenericMessenger* fMessenger;
		// The strength of the field
		G4double fB;
};

#endif
