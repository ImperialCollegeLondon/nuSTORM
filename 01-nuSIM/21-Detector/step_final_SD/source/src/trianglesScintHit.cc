// trianglesScintHit.cc
// Implementation of the trianglesScintHit class
//
#include "trianglesScintHit.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<trianglesScintHit>* trianglesScintHitAllocator=0;

// The constructor simply initializes all the values involved. fScintNb, fPlaneNb and fModuleNb
// are initialized as -1 because the copy numbers of logical volumes always start from 0.
trianglesScintHit::trianglesScintHit() : G4VHit(), fTrackID(-1), fScintNb(-1), fPlaneNb(-1), fModuleNb(-1), fEdep(0.0), fPos(G4ThreeVector()){}

// The destructor
trianglesScintHit::~trianglesScintHit(){}

// The other version of the constructior also initializes all the values involved, but using the value of right (whatever that is)
trianglesScintHit::trianglesScintHit(const trianglesScintHit &right) : G4VHit(), fTrackID(right.fTrackID), fScintNb(right.fScintNb), fEdep(right.fEdep), fPos(right.fPos){}

// Allows substitutions (the = operator) between two trianglesScintHit objects
const trianglesScintHit& trianglesScintHit::operator=(const trianglesScintHit &right){
	fTrackID = right.fTrackID;
	fScintNb = right.fScintNb;
	fPlaneNb = right.fPlaneNb;
	fModuleNb = right.fModuleNb;
	fEdep = right.fEdep;
	fPos = right.fPos;
	return *this;
}

G4bool trianglesScintHit::operator==(const trianglesScintHit &right) const{
	return (this == &right) ? true :false;
}

// Draw the hits (I think)
void trianglesScintHit::Draw(){
	auto visManager = G4VVisManager::GetConcreteInstance();
	if (! visManager) return;

	G4Circle circle(fPos);
	circle.SetScreenSize(4.0);
	circle.SetFillStyle(G4Circle::filled);
	G4Colour colour(1.0, 0.0, 0.0);
	G4VisAttributes attribs(colour);
	visManager->Draw(circle);
}
 void trianglesScintHit::Print(){
	 G4cout << "TrackID:" << fTrackID << "Strip number( Module, Plane, Rhombus, Scint):" << fModuleNb << fPlaneNb << fRhombNb << fScintNb << " Edep:" << std::setw(7) << G4BestUnit(fEdep, "Energy") << "Position:" <<std::setw(7)  <<G4BestUnit(fPos, "Length") <<G4endl;
 }
