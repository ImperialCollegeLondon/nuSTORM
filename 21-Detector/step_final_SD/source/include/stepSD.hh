// stepSD.hh
// Definition of the stepSD class.
// This class defines the sensitive detector that records position, energy and momentum for every step.
//
#ifndef stepSD_h
#define stepSD_h 1

#include "G4VSensitiveDetector.hh"

#include "stepHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

class stepSD : public G4VSensitiveDetector{
	public:
		stepSD(const G4String& name, const G4String& hitsCollectionName);
		virtual ~stepSD();

		virtual void Initialize(G4HCofThisEvent* hitsCollection);
		virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
		virtual void EndOfEvent(G4HCofThisEvent* hitsCollection);

	private:
		stepHitsCollection* fHitsCollection;
};

#endif
