// trianglesScintHit.hh
// Definition of the trianglesScintHit class.
// This class decides which quantities to record on a hit and provides the necessary variables.
// This class basically just defines the variables/containers, and trianglesScintSD class does the part where you actually retrieve the needed info from the hits.
//
#ifndef trianglesScintHit_h
#define trianglesScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"
//#include "G4LogicalVolume.hh"
//#include "G4Trancform3D.hh"
//#include "G4RotationMatrix.hh"

class G4ArrDef;
class G4AttValue;

// For now we try to record:
// - the strip number
// - the total energy deposit
//
//
class trianglesScintHit : public G4VHit{
	public:
		trianglesScintHit();
		trianglesScintHit(const trianglesScintHit&); // Two constructors depending on the input variables.
		virtual ~trianglesScintHit();

		const trianglesScintHit& operator=(const trianglesScintHit &right);  // These four lines allow you to use the = and == operators between trianglesScintHit objects
		G4bool operator==(const trianglesScintHit &right) const;             // in a similar way you would use with normal numbers, etc.
										     // They are not essential for the actual physics and the simulation.
		inline void *operator new(size_t);
		inline void operator delete(void*aHit);

		virtual void Draw();
		//virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
		//virtual std::vector<G4AttValue>* CreateAttValues() const;	(These two lines from example B5)
		virtual void Print();

		 // Method to handle data. Used for accumulating the total energy deposited to the scintillator planes.
		void Add(G4double de){
			fEdep += de;
		}

		// Methods to set the physical quantities in the instance 
		// In the entire application, we will get the values from the Step object and feed it to this to store it in the Hits collection.
		void SetTrackID (G4int track) {fTrackID = track;};
		void SetScintNb (G4int scint) {fScintNb = scint;};
		void SetRhombNb (G4int rhomb) {fRhombNb = rhomb;};
		void SetPlaneNb (G4int plane) {fPlaneNb = plane;};
		void SetModuleNb(G4int module){fModuleNb=module;};
		void SetEdep    (G4double de) { fEdep = de;};
		void SetPos     (G4ThreeVector xyz) {fPos = xyz;};

		// Methods to get the desired physical quantities from the instance.
		G4int GetTrackID() const {return fTrackID;};
		G4int GetScintNb() const {return fScintNb;};
		G4int GetRhombNb() const {return fRhombNb;};
		G4int GetPlaneNb() const {return fPlaneNb;};
		G4int GetModuleNb() const {return fModuleNb;};
		G4double GetEdep() const {return fEdep;};
		G4ThreeVector GetPos() const {return fPos;};

		private:
		// The values we try to store in a single Hit.
		G4int fTrackID; // The track ID(?)
		G4int fScintNb; // The Scintillator number.
		G4int fRhombNb; // The number of the rhombus that the scintillator belongs in.
		G4int fPlaneNb; // The number of the plane.
		G4int fModuleNb;// The number of the module.
		G4double fEdep; // The energy deposited.
		G4ThreeVector fPos; // The position of the event (Which coordinate?)
		
};

typedef G4THitsCollection<trianglesScintHit> trianglesScintHitsCollection;
// Use the G4HitsCollection<trianglesScintHit> class(?) to pack the Hits information into a HitsCollection.
// In terms of syntax, it's similar to when std::vector<MyType> creates a vector whose constituents are all MyType objects.

// Just copied the lines below from example B5/B2a.
extern G4ThreadLocal G4Allocator<trianglesScintHit>* trianglesScintHitAllocator;

inline void* trianglesScintHit::operator new(size_t){
	if (!trianglesScintHitAllocator){
		trianglesScintHitAllocator = new G4Allocator<trianglesScintHit>;
	}
	return (void*)trianglesScintHitAllocator->MallocSingle();
}

inline void trianglesScintHit::operator delete(void*aHit){
	trianglesScintHitAllocator->FreeSingle((trianglesScintHit*) aHit);
}

#endif
