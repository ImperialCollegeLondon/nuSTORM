//trianglesPrimaryGeneratorAction.cc
//Implementation of th trianglesPrimaryGeneratorAction class. Definitions should be found in (project source)/src/trianglesPrimaryGeneratorAction.hh
//
#include "trianglesPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

// The constructor of the class. Generates the particle gun as the instance is created.
trianglesPrimaryGeneratorAction::trianglesPrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(0){
	// Generates a particle gun that shoots one particle at a time (I think).
	G4int n_particle = 1;
	fParticleGun = new G4ParticleGun(n_particle);

	// The kinematic properties of the particle being shot.
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle = particleTable->FindParticle(particleName="mu+");
	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
	fParticleGun->SetParticleEnergy(0.5*GeV);
}

// The destructor. Deletes (frees) the particle gun produced by the constructor.
trianglesPrimaryGeneratorAction::~trianglesPrimaryGeneratorAction(){
	delete fParticleGun;
}

void trianglesPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
	G4cout << "trianglesPrimaryGeneratorAction::GeneratePrimaries" << G4endl;
	G4double world_halfXY = 0;
	G4double world_halfZ = 0;

	auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
	G4Box* worldBox = nullptr;
	if (worldLV){
		worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
	}
	if (worldBox){
		world_halfXY =0.4 * worldBox->GetXHalfLength();
		world_halfZ = worldBox->GetZHalfLength();
	}
	else {
		G4ExceptionDescription msg;
		msg << "World volume of box shape not found." << G4endl;
		msg << "The gun will be placed in the center.";
		G4Exception("trianglesPrimaryGeneratorAction::GeneratePrimaries()", "MyCode0002", JustWarning, msg);
	}

	fParticleGun->SetParticlePosition(G4ThreeVector(world_halfXY*(2*G4UniformRand()-1),world_halfXY*(2*G4UniformRand()-1),world_halfZ*(2*G4UniformRand()-1)));
	fParticleGun->SetParticleEnergy((0.5+4.5*G4UniformRand())*GeV);
	fParticleGun->GeneratePrimaryVertex(anEvent);
}
