// stepSD.cc
// Implementation of the stepSD class.
//
#include "stepHit.hh"
#include "stepSD.hh"
#include "trianglesConstants.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"


stepSD::stepSD(const G4String& name, const G4String& hitsCollectionName) : G4VSensitiveDetector(name), fHitsCollection(NULL){
	collectionName.insert(hitsCollectionName);
}

stepSD::~stepSD() {}

void stepSD::Initialize(G4HCofThisEvent* hitsCollection){
	G4cout << "stepSD::Initialize" << G4endl;
	fHitsCollection = new stepHitsCollection(SensitiveDetectorName, collectionName[0]);
	
	G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

	hitsCollection->AddHitsCollection(hcID, fHitsCollection);
}

G4bool stepSD::ProcessHits(G4Step* aStep, G4TouchableHistory*){
	auto edep = aStep->GetTotalEnergyDeposit();
	//G4cout << "stepSD::ProcessHits, edep = " << edep << G4endl;
	if (edep==0.0) return true;

	// Get the strip in which the step took place.
	auto preStepPoint = aStep->GetPreStepPoint();
	auto touchable = preStepPoint->GetTouchable();
	auto scintNo = touchable->GetVolume()->GetCopyNo();
	auto rhombNo = touchable->GetVolume(1)->GetCopyNo();
	auto planeNo = touchable->GetVolume(2)->GetCopyNo();
	auto moduleNo = touchable->GetVolume(3)->GetCopyNo();
	//auto stripNo = moduleNo*6*kNofScintillators + planeNo*kNofScintillators + rhombNo;

	// Get the position and time of the pre-step point.
	auto hitTime = preStepPoint->GetGlobalTime();
	auto worldPos = preStepPoint->GetPosition();
	auto momentumDir = preStepPoint->GetMomentumDirection();
	auto momentumVec = preStepPoint->GetMomentum();

	// Generate a hit.
	auto hit = new stepHit();
	hit->SetScintNb(scintNo);
	hit->SetRhombNb(rhombNo);
	hit->SetPlaneNb(planeNo);
	hit->SetModuleNb(moduleNo);
	hit->SetGlobalTime(hitTime);
	hit->SetEdep(edep);
	hit->SetPos(worldPos);
	hit->SetMomentumDir(momentumDir);
	hit->SetMomentumMag(momentumVec.mag());
	hit->SetMomentumVec(momentumVec);

	fHitsCollection->insert(hit);

	return true;
}

void stepSD::EndOfEvent(G4HCofThisEvent*){
	G4cout << "stepSD::EndOfEvent" << G4endl;
	if (verboseLevel>1){
		G4int nofHits = fHitsCollection->entries();
		G4cout << G4endl
			<<"------>Hits Collection: in this event there were " << nofHits
			<< "Hits (step) in the scintillator planes:" << G4endl;
		for (G4int i=0; i<nofHits; i++) (*fHitsCollection)[i]->Print();
	}
}
