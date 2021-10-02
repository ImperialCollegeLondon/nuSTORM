// trianglesPrimaryGeneratorAction.hh
// Brief definiton of the trianglesPrimaryGeneratorAction class.
// This class defines where and how the original particle beams (or the 'events') are generated before entering the detector (of course it can be inside the detector as well).
//
// This class is required for the implementation of the trianglesActionInitialization class.
//
#ifndef trianglesPrimaryGeneratorAction_h
#define trianglesPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class G4Box;

class trianglesPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction{

	public:
		trianglesPrimaryGeneratorAction();
		virtual ~trianglesPrimaryGeneratorAction();

		// method from the base class G4VUserPrimaryGeneratorAction. Inplement this on another file to define the details/initialization of the events (e.g. type and energy of particle)
		virtual void GeneratePrimaries(G4Event*);

		//method to access the particle gun
		const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

	private:
		G4ParticleGun* fParticleGun; //pointer to the G4 gun class
		//G4Box* fEnvelopeBox;
};

#endif
