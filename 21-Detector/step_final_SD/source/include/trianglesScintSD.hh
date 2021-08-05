// trianglesScintSD.hh
// Definition of the trianglesScintSD class.
// The Hits generated according to the trianglesScintHit class are accounted in this class.
// This is where we actually take out the information we want from the hits.
//
//
#ifndef trianglesScintSD_h
#define trianglesScintSD_h 1

#include "G4VSensitiveDetector.hh"

#include "trianglesScintHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent; // Geant4 Hit Collection of this event

class trianglesScintSD : public G4VSensitiveDetector{
	public:
		trianglesScintSD(const G4String& name, const G4String& hitsCollectionName, G4int num_strips);
		// Declare the name and number of scintillator strips when declaring this sensitive detector.
		virtual ~trianglesScintSD();

		// Method from the base class, i.e. the G4VSensitiveDetector class
		// Initialize the hits collection at the beginning of an event(I think)
		virtual void Initialize(G4HCofThisEvent* hitCollection);
		// Process the Hits at each Step.
		virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
		// Sum up all the hits at the end of the event
		virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

	private:
		trianglesScintHitsCollection* fHitsCollection;
		// The collection of all the hits.
		G4int fNofStrips;
		// The number of scintillator strips.
};

#endif
