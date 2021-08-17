// trianglesRunAction.cc
// Implementation of the trianglesRunAction class.
//
#include "trianglesRunAction.hh"
#include "trianglesEventAction.hh"
#include "trianglesConstants.hh"

#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4UIcommand.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericAnalysisManager.hh"
//#include "g4root.hh"

using G4AnalysisManager = G4GenericAnalysisManager;

// The constructor.
trianglesRunAction::trianglesRunAction(trianglesEventAction* eventAction) : G4UserRunAction(), fEventAction(eventAction){
	// Create the generic analysis manager
	// The choice of analysis technology (e.g. root, csv, ...) is done according to the file extension.
	// If no extension is given, Root output would be chosen as the default.
	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetVerboseLevel(1);

	//Default settings
	analysisManager->SetNtupleMerging(true);
	// This allows the merging of ntuples (available only with Root output)
	analysisManager->SetFileName("triangleScints");

	// Creating 1D Histograms
	analysisManager->CreateH1("wholeEnergy", "Energy Deposited", 50, 0., 1000); // (name, title, number of bins, minimum, maximum)

	// Create an ntuple
	if (fEventAction){
		analysisManager->CreateNtuple("trianglesScints", "Energy");
		analysisManager->CreateNtupleIColumn("EventNo"); // Column ID 0
		analysisManager->CreateNtupleIColumn("InitParticleType"); // ID 1 
		analysisManager->CreateNtupleDColumn("InitialEnergy"); //ID 2
		analysisManager->CreateNtupleDColumn("InitialMomX"); //ID 3
		analysisManager->CreateNtupleDColumn("InitialMomY"); //ID 4
		analysisManager->CreateNtupleDColumn("InitialMomZ"); //ID 5
		analysisManager->CreateNtupleDColumn("InitialPosX"); //ID 6
		analysisManager->CreateNtupleDColumn("InitialPosY"); // ID 7
		analysisManager->CreateNtupleDColumn("InitialPosZ"); // ID 8
		analysisManager->CreateNtupleIColumn("StripNo_ET", eventAction->GetStripNo()); // ID 9
		analysisManager->CreateNtupleIColumn("PlaneNo_ET", eventAction->GetPlaneNo()); // ID 10
		analysisManager->CreateNtupleIColumn("GlobalPlaneNo_ET", eventAction->GetPlaneNoGlobal()); // ID 11
		analysisManager->CreateNtupleIColumn("ModuleNo_ET", eventAction->GetModuleNo()); // ID 12
		analysisManager->CreateNtupleDColumn("EnergyDeposit_ET", eventAction->GetEdep()); // ID 13
		analysisManager->CreateNtupleDColumn("TotalEnergy_ET"); // ID 14

		analysisManager->CreateNtupleIColumn("StripNo_step", eventAction->GetSStripNo()); // ID 15
		analysisManager->CreateNtupleIColumn("PlaneNo_step", eventAction->GetSPlaneNo()); // ID 16
		analysisManager->CreateNtupleIColumn("ModuleNo_step", eventAction->GetSModuleNo()); // ID 17
		analysisManager->CreateNtupleDColumn("EDep_step", eventAction->GetSEdep()); // ID 18
		analysisManager->CreateNtupleDColumn("globalPosX_step", eventAction->GetSPosX()); // ID 19
		analysisManager->CreateNtupleDColumn("globalPosY_step", eventAction->GetSPosY()); // ID 20
		analysisManager->CreateNtupleDColumn("globalPosZ_step", eventAction->GetSPosZ()); // ID 21
		analysisManager->CreateNtupleDColumn("momMag_step", eventAction->GetSMomMag()); // ID 22
		analysisManager->CreateNtupleDColumn("momX_step", eventAction->GetSMomX()); // ID 23
		analysisManager->CreateNtupleDColumn("momY_step", eventAction->GetSMomY()); // ID 24
		analysisManager->CreateNtupleDColumn("momZ_step", eventAction->GetSMomZ()); // ID 25
		analysisManager->CreateNtupleDColumn("time_step", eventAction->GetSGlobalTime()); // ID 26
		// Let's think about the momentum direction later.
		analysisManager->FinishNtuple();
		//ET stands for event total. For quantities that are summed up throughout the event.
	}

	analysisManager->SetNtupleFileName(0, "Scints_ntuple");
}

trianglesRunAction::~trianglesRunAction(){
	delete G4AnalysisManager::Instance();
}

void trianglesRunAction::BeginOfRunAction(const G4Run*){

	// Get the analysis manager
	auto analysisManager = G4AnalysisManager::Instance();

	analysisManager->OpenFile();
}

void trianglesRunAction::EndOfRunAction(const G4Run*){
	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->Write();
	analysisManager->CloseFile();
}

