// trianglesEventAction.hh
// Definition of the trianglesEventAction class.
// Defines the actions that will be taken with each of the events.
//
#ifndef trianglesEventAction_h
#define trianglesEventAction_h 1

#include "G4UserEventAction.hh"

#include "trianglesScintHit.hh"
#include "stepHit.hh"

#include "globals.hh"


class trianglesEventAction : public G4UserEventAction{
	public:
		trianglesEventAction();
		virtual ~trianglesEventAction();

		virtual void BeginOfEventAction(const G4Event*);
		virtual void EndOfEventAction(const G4Event*);

		std::vector<G4int>& GetStripNo(){return fStripNoVec;}
		std::vector<G4int>& GetPlaneNo(){return fPlaneNoVec;}
		std::vector<G4int>& GetPlaneNoGlobal(){return fPlaneNoGlobalVec;}
		std::vector<G4int>& GetModuleNo(){return fModuleNoVec;}
		std::vector<G4double>& GetEdep(){return fEdepVec;}

		std::vector<G4int>& GetSStripNo(){return fS_StripNoVec;};
		std::vector<G4int>& GetSPlaneNo(){return fS_PlaneNoVec;};
		std::vector<G4int>& GetSPlaneNoGlobal(){return fS_PlaneNoGlobalVec;};
		std::vector<G4int>& GetSModuleNo(){return fS_ModuleNoVec;};
		std::vector<G4double>& GetSEdep(){return fS_EdepVec;};
		std::vector<G4double>& GetSPosX(){return fS_posX;};
		std::vector<G4double>& GetSPosY(){return fS_posY;};
		std::vector<G4double>& GetSPosZ(){return fS_posZ;};
		std::vector<G4ThreeVector>& GetSMomDit() {return fS_momDir;};
		std::vector<G4double>& GetSMomMag() {return fS_momMag;};
		std::vector<G4double>& GetSMomX() {return fS_momX;};
		std::vector<G4double>& GetSMomY() {return fS_momY;};
		std::vector<G4double>& GetSMomZ() {return fS_momZ;};
		std::vector<G4double>& GetSGlobalTime(){return fS_time;};

	private:
		trianglesScintHitsCollection* GetHitsCollection(G4int hcID, const G4Event* event) const; 
		stepHitsCollection* GetSHitsCollection(G4int hcID, const G4Event* event) const;
		//void PrintEventStatistics(G4double scintEdep) const;

		G4int fScintHCID; // The ID of the Scintillator Hits collection.
		G4int fScint2HCID; // ID of the Scintillator2 Hits collection.
		G4int fStepHCID;  // ID of the step (sensitive detector) Hits collection.
		G4int fStep2HCID; // ID of the step2 Hits collection.

		std::vector<G4int> fStripNoVec; // Vector to store the scintillator strip no that fired.
		std::vector<G4int> fPlaneNoVec; // Vector to store the detector plane number (the one that ranges from 0 to 5).
		std::vector<G4int> fPlaneNoGlobalVec; // Vector to store the global detector plane number;
		std::vector<G4int> fModuleNoVec; // Vector to store the module number.
		std::vector<G4double> fEdepVec; // Vector to store the energy deposit.

		std::vector<G4int> fS_StripNoVec; // Vector to store the strip no. step by step.
		std::vector<G4int> fS_PlaneNoVec;
		std::vector<G4int> fS_PlaneNoGlobalVec;
		std::vector<G4int> fS_ModuleNoVec;
		std::vector<G4double> fS_EdepVec;
		std::vector<G4double> fS_posX;
		std::vector<G4double> fS_posY;
		std::vector<G4double> fS_posZ;
		std::vector<G4ThreeVector> fS_momDir;
		std::vector<G4double> fS_momMag;
		std::vector<G4double> fS_momX;
		std::vector<G4double> fS_momY;
		std::vector<G4double> fS_momZ;
		std::vector<G4double> fS_time;
};

#endif
