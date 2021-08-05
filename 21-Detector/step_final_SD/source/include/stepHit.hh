
// stepHit.hh
// Definition of the stepHit class.
// This is a class to define a new sensitive detector that records the energy, momentum and position for each step.
//
#ifndef stepHit_h
#define stepHit_h 1

#include "trianglesConstants.hh"

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

class G4ArrDef;
class G4AttValue;

class stepHit : public G4VHit{
	public:
		stepHit();
		stepHit(const stepHit&);
		virtual ~stepHit();

		const stepHit& operator=(const stepHit &right);
		G4bool operator==(const stepHit &right) const;

		inline void *operator new(size_t);
		inline void operator delete(void*aHit);

		//virtual void Draw(); This would prabably make the canvas a bit too crowded.
		virtual void Print();

		void SetTrackID (G4int track) {fTrackID = track;};
		void SetScintNb (G4int scint) {fScintNb = scint;};
		void SetRhombNb (G4int rhomb) {fRhombNb = rhomb;};
		void SetPlaneNb (G4int plane) {fPlaneNb = plane;};
		void SetModuleNb (G4int module) {fModuleNb = module;};
		void SetGlobalTime(G4double time) {fGlobalTime = time;};
		void SetEdep (G4double edep) {fEdep = edep;};
		void SetPos (G4ThreeVector pos) {fPos = pos;};
		void SetMomentumDir (G4ThreeVector momentum) {fMomentumDir = momentum;};
		void SetMomentumMag (G4double momentumMag) {fMomentumMag = momentumMag;};
		void SetMomentumVec (G4ThreeVector momentum) {fMomentumVec = momentum;};

		G4int GetTrackID() const {return fTrackID;};
		G4int GetScintNb () const {return fScintNb;};
		G4int GetRhombNb () const {return fRhombNb;};
		G4int GetPlaneNb () const {return fPlaneNb;};
		G4int GetModuleNb () const {return fModuleNb;};
		G4double GetGlobalTime() const {return fGlobalTime;};
		G4double GetEdep() const {return fEdep;};
		G4ThreeVector GetPos() const {return fPos;};
		G4ThreeVector GetMomentumDir() const {return fMomentumDir;};
		G4double GetMomentumMag() const {return fMomentumMag;};
		G4ThreeVector GetMomentumVec() const {return fMomentumVec;};
	
		//G4int GetStripID() const{return fModuleNb*6*kNofScintillators+fPlaneNb*kNofScintillators+fRhombNb;};
	private:
		G4int fTrackID;
		G4int fScintNb;
		G4int fRhombNb;
		G4int fPlaneNb;
		G4int fModuleNb;
		G4double fEdep;
		G4double fGlobalTime; // Time from the beginning of the event to the pre-step point.
		G4ThreeVector fPos; // Global position at the pre-step point (i.e. position relative to the world)
		G4ThreeVector fMomentumDir; // Direction of the momentum at the pre-step point (unit vector)
		G4double fMomentumMag; // Magnitude of the momentum at the pre-step point
		G4ThreeVector fMomentumVec; // Momentum vector
};

typedef G4THitsCollection<stepHit> stepHitsCollection;

extern G4ThreadLocal G4Allocator<stepHit>* stepHitAllocator;

inline void* stepHit::operator new(size_t){
	if (!stepHitAllocator){
		stepHitAllocator = new G4Allocator<stepHit>;
	}
	return (void*)stepHitAllocator->MallocSingle();
}

inline void stepHit::operator delete(void*aHit){
	stepHitAllocator->FreeSingle((stepHit*) aHit);
}
#endif
