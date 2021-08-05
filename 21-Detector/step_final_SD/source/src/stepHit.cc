// stepHit.cc
// IMplementation of the class stepHit
//
#include "stepHit.hh"

#include "G4UnitsTable.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<stepHit>* stepHitAllocator=0;

stepHit::stepHit() : G4VHit(),
	fTrackID(-1),
	fScintNb(-1),
	fRhombNb(-1),
	fPlaneNb(-1),
	fModuleNb(-1),
	fEdep(0.0),
	fGlobalTime(0.0),
	fPos(G4ThreeVector()),
	fMomentumDir(G4ThreeVector()),
	fMomentumMag(0.0),
	fMomentumVec(G4ThreeVector())
{}

stepHit::~stepHit(){}

stepHit::stepHit(const stepHit &right) : G4VHit(),
	fTrackID(right.fTrackID),
	fScintNb(right.fScintNb),
	fRhombNb(right.fRhombNb),
	fPlaneNb(right.fPlaneNb),
	fModuleNb(right.fModuleNb),
	fEdep(right.fEdep),
	fGlobalTime(right.fGlobalTime),
	fPos(right.fPos),
	fMomentumDir(right.fMomentumDir),
	fMomentumMag(right.fMomentumMag),
	fMomentumVec(right.fMomentumVec)
{}

const stepHit& stepHit::operator=(const stepHit &right){
	fTrackID = right.fTrackID;
	fScintNb = right.fScintNb;
	fRhombNb = right.fRhombNb;
	fPlaneNb = right.fPlaneNb;
	fModuleNb = right.fModuleNb;
	fEdep = right.fEdep;
	fGlobalTime = right.fGlobalTime;
	fPos = right.fPos;
	fMomentumDir = right.fMomentumDir;
	fMomentumMag = right.fMomentumMag;
	fMomentumVec = right.fMomentumVec;
	return *this;
}

G4bool stepHit::operator==(const stepHit &right) const{
	return (this == &right) ? true : false;
}

void stepHit::Print(){
	G4cout << "TrackID:" << fTrackID << "Strip Number:" << fModuleNb << fPlaneNb << fRhombNb << fScintNb << "Edep:" << std::setw(7) << G4BestUnit(fEdep, "Energy") << "Position:" << std::setw(7) << G4BestUnit(fPos, "Length") << G4endl;
}
