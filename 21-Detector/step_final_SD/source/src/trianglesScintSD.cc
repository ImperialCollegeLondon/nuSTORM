// trianglesScintSD.cc
// Implementation of the trianlgesScintSd class.
//
//

#include "trianglesScintSD.hh"
#include "trianglesScintHit.hh"
#include "trianglesConstants.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"


trianglesScintSD::trianglesScintSD(const G4String& name, const G4String& hitsCollectionName, G4int num_strips) : G4VSensitiveDetector(name), fHitsCollection(NULL), fNofStrips(num_strips){
	// Define the name of the hits collection that this sensitive detetctor will handle.
	// The G4String& name is the name used when declaring an instance
	// In our case, this is done in trianglesDetectorConstruction::createSDandField.

	// nameColl = /scintillatorsColl, /scintillators2Coll
	collectionName.insert(hitsCollectionName);

       	// collectionName is a vector defined in the base class G4VSensitiveDetector. It must be declared(?) in the constructor
	// If a sensitive detector produces multiple number of hits, we have to do collectionName.insert(name) multiple times. 
	// In that case, the variable collectionName will become a vector with multiple components. In our case, it only has the 0th complnent.
}

trianglesScintSD::~trianglesScintSD() {}

void trianglesScintSD::Initialize(G4HCofThisEvent* hce){
	G4cout << "trianglesScintSD::Initialize" << G4endl;
	// This action is called at the beginning of every event.
	// Create the hits collection (fHitsCollection is a member variable of this (trianglesScintSD) class)
	fHitsCollection = new trianglesScintHitsCollection(SensitiveDetectorName, collectionName[0]);

	// We register the hits collection we just made to the G4SDManager. First we get the collection ID...
	G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	// and then we add it to the G4HCofThisEvent object.
	hce->AddHitsCollection(hcID, fHitsCollection); 
       
	// Now the hits collection has a name, name of the logical volume, and a hits collection ID registered.
	//
	// Next, we generate a hit corresponding to each of the scintillator strips.
	for (G4int i=0; i<fNofStrips+1; i++){
		fHitsCollection->insert(new trianglesScintHit());
	}
	G4cout << "scint Hits collection created, entries = " << fHitsCollection->entries() << G4endl;
}

G4bool trianglesScintSD::ProcessHits(G4Step* aStep, G4TouchableHistory*){
	// This function first receives a step. Then, it will try to locate the scintillator strip that the step took place in (hitID).
	// After that, it will choose the appropriate hit in the hits collection(*fHitsCollection)[hitID], and add the energy deposited in the strip, and the entire detector.
	// After a bunch of steps, the nth component of fHitsCollection would have the total energy the nth scintillator strip has received,
	// and the last component would have the total energy the entire scintillator detector has received.
	
	// The energy deposited in this step.
	G4double edep = aStep->GetTotalEnergyDeposit();
	//G4cout << "trianglesScintSD::ProcessHits, edep = " << edep << G4endl;

	if (edep==0.0) return false;

	// Get the hits in the scintillator strips.
	auto touchable = aStep->GetPreStepPoint()->GetTouchable();
	//auto scintNb = touchable->GetCopyNumber();
	auto rhombNb = touchable->GetCopyNumber(1);
	auto planeNb = touchable->GetCopyNumber(2);
	auto moduleNb= touchable->GetCopyNumber(3);
	auto hitID = rhombNb+kNofScintillators*planeNb+moduleNb*6*kNofScintillators; 
	auto hit = (*fHitsCollection)[hitID];
	if (!hit){
		G4ExceptionDescription msg;
		msg << "Cannot access hit "<< hitID;
		G4Exception("trianglesScintSD::ProcessHits()", "MyCode0004", FatalException, msg);
	}

	// get the hit for accounting of the total energy
	auto hitTotal = (*fHitsCollection)[fHitsCollection->entries()-1];

	hit->Add(edep);
	hitTotal->Add(edep);


	return true;
}

void trianglesScintSD::EndOfEvent(G4HCofThisEvent*){
	G4cout << "trianglesScintSD::EndOfEvent" << G4endl;
	if (verboseLevel>1) {
		G4int nofHits = fHitsCollection->entries();
		G4cout<< G4endl
		      <<"-------->Hits Collection: in this event there were " << nofHits
		      <<" hits in the scintillator planes: " <<G4endl;
		for (G4int i=0; i<nofHits; i++) (*fHitsCollection)[i]->Print();
	}
}
