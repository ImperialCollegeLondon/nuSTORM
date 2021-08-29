//trianglesEventAction.cc
//
//Implementation of the trianglesEventAction class.
//
#include "trianglesEventAction.hh"
#include "trianglesScintHit.hh"
#include "trianglesConstants.hh"
#include "stepHit.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "g4analysis.hh"

using std::array;
using std::vector;

/* This was part of the example B5 code, but I don't think we need it for our purpose because we only have one sensitive detector (i.e. only one set of hit colletion).
namespace {
	// Function GetHC takes an event and a collID, and returns a hit collection with the given ID.
	G4VHitsCollection* GetHC(const G4Event* event, G4int collId){
		auto hce = event->GetHCofThisEvent();
		if (!hce) {
			G4ExceptionDescription msg;
			msg << "No hits collection of this event found." << G4endl;
			G4Exception("trianglesEventAction::EndOfEventAction()", "trianglesCode001", JustWarning, msg);
			return nullptr;
		}

		auto hc = hce->GetHC(collId);
		if (!hc){
			G4ExceptionDescription msg;
			msg << "HitsCollection" << collId << "of this event not found." << G4endl;
			G4Exception("trianglesEventAction::EndOfEventAction()", "trianglesCode001", JustWarning, msg);
		}
		return hc;
	}
}

*/

trianglesEventAction::trianglesEventAction() :G4UserEventAction(), 
	fScintHCID(-1), 
	fScint2HCID(-1),
	fStepHCID(-1),
	fStep2HCID(-1) ,
	fStripNoVec(0), 
	fPlaneNoVec(0), 
	fPlaneNoGlobalVec(0), 
	fModuleNoVec(0), 
	fEdepVec(0), 
	fS_StripNoVec(0),
	fS_PlaneNoVec(0),
	fS_PlaneNoGlobalVec(0),
	fS_ModuleNoVec(0),
	fS_EdepVec(0),
	fS_posX(0),
	fS_posY(0),
	fS_posZ(0),
	fS_momDir(0),
	fS_momMag(0),
	fS_time(0)
{
	G4RunManager::GetRunManager()->SetPrintProgress(1);
}

trianglesEventAction::~trianglesEventAction() {}

trianglesScintHitsCollection* trianglesEventAction::GetHitsCollection(G4int hcID, const G4Event* event) const {
	auto hitsCollection = static_cast<trianglesScintHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

	if (!hitsCollection){
		G4ExceptionDescription msg;
		msg << "Cannot access hitsCollection ID" << hcID;
		G4Exception("trianglesEventAction::GetHitsCollection()", "MyCode0003", FatalException, msg);
	}
	return hitsCollection;
}

stepHitsCollection* trianglesEventAction::GetSHitsCollection(G4int hcID, const G4Event*event) const{
	auto hitsCollection = static_cast<stepHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

	if (!hitsCollection){
		G4ExceptionDescription msg;
		msg << "Cannot access Step hitsCollection ID" << hcID;
		G4Exception("trianglesEventAction::GetSHitsCollection()", "MyCode0003", FatalException, msg);
	}
	return hitsCollection;
}

void trianglesEventAction::BeginOfEventAction(const G4Event*) {
	G4cout << "trianglesEventAction::BeginOfEventAction" << G4endl;
	fStripNoVec = {};
	fPlaneNoVec = {};
	fPlaneNoGlobalVec = {};
	fModuleNoVec = {};
	fEdepVec = {};
	fS_StripNoVec = {};
	fS_PlaneNoVec = {};
	fS_PlaneNoGlobalVec = {};
	fS_ModuleNoVec = {};
	fS_EdepVec = {};
	fS_posX = {};
	fS_posY = {};
	fS_posZ = {};
	fS_momDir = {};
	fS_momMag = {};
	fS_time = {};

}

void trianglesEventAction::EndOfEventAction(const G4Event* event){
	G4cout << "trianglesEventAction::EndOfEventAction" << G4endl;
	// First, access the event and retrieve a hits collection once.
	if (fScintHCID == -1) {
		G4cout << "Getting Hits Collection IDs" << G4endl;
		fScintHCID = G4SDManager::GetSDMpointer()->GetCollectionID("scintColl");
		fScint2HCID= G4SDManager::GetSDMpointer()->GetCollectionID("scint2Coll");
		fStepHCID = G4SDManager::GetSDMpointer()->GetCollectionID("scint_stepsColl");
		fStep2HCID = G4SDManager::GetSDMpointer()->GetCollectionID("scint2_stepsColl");
	       	// Get the hits collection with the names 'scintColl' and 'scint2Coll'
		// The names were defined when the trianglesScintSD object was created in the trianglesDetectorConstruction class.
		G4cout << "fScintHCID=" << fScintHCID << ", fScint2HCID=" << fScint2HCID << ", fStepHCID=" << fStepHCID << ", fStep2HCID=" << fStep2HCID << G4endl;
	}
	// Next, get the hits collections according to the ID that we just got.
	auto scintHC = GetHitsCollection(fScintHCID, event);
	auto scint2HC = GetHitsCollection(fScint2HCID, event);
	G4int eventID = event->GetEventID();

	G4int primaryPartID = event->GetPrimaryVertex(0)->GetPrimary(0)->GetPDGcode();
	auto primaryPos = event->GetPrimaryVertex(0)->GetPosition();

	auto primaryEnergy = event->GetPrimaryVertex(0)->GetPrimary(0)->GetKineticEnergy();
	auto primaryMom = event->GetPrimaryVertex(0)->GetPrimary(0)->GetMomentum();
	/*
	// Print the event summary (probably needs some modifications)
	//auto eventID = event->GetEventID();
	auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
	if ((printModulo > 0) && (eventID % printModulo == 0)){
		G4cout<<"---> End of Event:" << eventID << G4endl;

		PrintEventStatistics(scintHC[scintHC->entries()]->GetEdep());
	}
	*/

	// Get the hits and fill the histograms & ntuple
	//
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	analysisManager->FillNtupleIColumn(0, eventID); // Fill the first column with the event ID.
	analysisManager->FillNtupleIColumn(1, primaryPartID);
	analysisManager->FillNtupleDColumn(2, primaryEnergy);
	analysisManager->FillNtupleDColumn(3, primaryMom.getX());
	analysisManager->FillNtupleDColumn(4, primaryMom.getY());
	analysisManager->FillNtupleDColumn(5, primaryMom.getZ());
	analysisManager->FillNtupleDColumn(6, primaryPos.getX());
	analysisManager->FillNtupleDColumn(7, primaryPos.getY());
	analysisManager->FillNtupleDColumn(8, primaryPos.getZ());

	// Get the hits corresponding to each of the scintillator strips.
	for (G4int j=0; j<kNofStrips; j++){
		auto scintHit = (*scintHC)[j];
		if( scintHit->GetEdep()){
			G4int modNo = j/(6*kNofScintillators);
			G4int planeNo = (j-modNo*6*kNofScintillators)/kNofScintillators;
			G4int stripNo = 2* (j-modNo*6*kNofScintillators-planeNo*kNofScintillators);
			G4int planeNoGlobal = 6*modNo + planeNo;

			fStripNoVec.emplace_back(stripNo);
			fPlaneNoVec.emplace_back(planeNo);
			fPlaneNoGlobalVec.emplace_back(planeNoGlobal);
			fModuleNoVec.emplace_back(modNo);
			fEdepVec.emplace_back(scintHit->GetEdep());
		}
	}
	for (G4int j=0; j<kNofStrips; j++){
		auto scint2Hit = (*scint2HC)[j];
		if( scint2Hit->GetEdep()){
			G4int modNo = j/(6*kNofScintillators);
			G4int planeNo = (j-modNo*6*kNofScintillators)/kNofScintillators;
			G4int stripNo = 2* (j-modNo*6*kNofScintillators-planeNo*kNofScintillators) + 1;
			G4int planeNoGlobal = 6*modNo + planeNo;

			fStripNoVec.emplace_back(stripNo);
			fPlaneNoVec.emplace_back(planeNo);
			fPlaneNoGlobalVec.emplace_back(planeNoGlobal);
			fModuleNoVec.emplace_back(modNo);
			fEdepVec.emplace_back(scint2Hit->GetEdep());
		}
	}
	// Get the final hit in the hits collection. This hit contains the total energy deposited to the entire detector.
	auto scint_entireHit = (*scintHC)[scintHC->entries()-1]; // If everything is coded correctly, scintHC->entries() = kNofStrips + 1
	auto scint2_entireHit= (*scint2HC)[scint2HC->entries()-1];

	G4double entireEdep = scint_entireHit->GetEdep() + scint2_entireHit->GetEdep();
	analysisManager->FillH1(0, entireEdep); // Total energy deposited in the scintillators.
	analysisManager->FillNtupleDColumn(14, entireEdep); // Total energy deposited in all of the scintillators.

	// Next, the part which records quantities per each step.
	auto StepHC = GetSHitsCollection(fStepHCID, event);
	auto Step2HC = GetSHitsCollection(fStep2HCID, event);

	//G4int stepNum = StepHC->entries();
	//G4int step2Num = Step2HC->entries();
	for (unsigned long i=0; i<StepHC->entries(); i++){
		auto stepHit = (*StepHC)[i];
		fS_StripNoVec.emplace_back(2*stepHit->GetRhombNb());
		fS_PlaneNoVec.emplace_back(stepHit->GetPlaneNb());
		fS_PlaneNoGlobalVec.emplace_back(stepHit->GetPlaneNb()+6*stepHit->GetModuleNb());
		fS_ModuleNoVec.emplace_back(stepHit->GetModuleNb());
		fS_EdepVec.emplace_back(stepHit->GetEdep());
		fS_posX.emplace_back(stepHit->GetPos().getX());
		fS_posY.emplace_back(stepHit->GetPos().getY());
		fS_posZ.emplace_back(stepHit->GetPos().getZ());
		fS_momMag.emplace_back(stepHit->GetMomentumMag());
		fS_momX.emplace_back(stepHit->GetMomentumVec().getX());
		fS_momY.emplace_back(stepHit->GetMomentumVec().getY());
		fS_momZ.emplace_back(stepHit->GetMomentumVec().getZ());
		fS_time.emplace_back(stepHit->GetGlobalTime());
	}
	for (unsigned long i=0; i<Step2HC->entries(); i++){
		auto step2Hit = (*Step2HC)[i];
		fS_StripNoVec.emplace_back(2*step2Hit->GetRhombNb()+1);
		fS_PlaneNoVec.emplace_back(step2Hit->GetPlaneNb());
		fS_PlaneNoGlobalVec.emplace_back(step2Hit->GetPlaneNb()+6*step2Hit->GetModuleNb());
		fS_ModuleNoVec.emplace_back(step2Hit->GetModuleNb());
		fS_EdepVec.emplace_back(step2Hit->GetEdep());
		fS_posX.emplace_back(step2Hit->GetPos().getX());
		fS_posY.emplace_back(step2Hit->GetPos().getY());
		fS_posZ.emplace_back(step2Hit->GetPos().getZ());
		fS_momMag.emplace_back(step2Hit->GetMomentumMag());
		fS_momX.emplace_back(step2Hit->GetMomentumVec().getX());
		fS_momY.emplace_back(step2Hit->GetMomentumVec().getY());
		fS_momZ.emplace_back(step2Hit->GetMomentumVec().getZ());
		fS_time.emplace_back(step2Hit->GetGlobalTime());
	}
	analysisManager->AddNtupleRow();
}

