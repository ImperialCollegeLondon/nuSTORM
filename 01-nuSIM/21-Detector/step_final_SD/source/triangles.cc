// triangles.cc
// The main function of the simulation.
//
#include "trianglesDetectorConstruction.hh"
#include "trianglesActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

int main(int argc, char** argv){

	//Detect interactive mode and define UI session
	//
	G4UIExecutive* ui = 0;
	if (argc == 1) {
		ui = new G4UIExecutive(argc, argv);
	}

	// Construct the default rum manager
	//
	auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

	//Set the mandatory initialization classes
	//
	//Detector construction
	runManager->SetUserInitialization(new trianglesDetectorConstruction);

	// List of the physics involved (i.e. the interactions between the incident particle and the detector volume)
	auto physlist = new FTFP_BERT;
	physlist->RegisterPhysics(new G4StepLimiterPhysics());
	runManager->SetUserInitialization(physlist);

	// User action initialization
	runManager->SetUserInitialization(new trianglesActionInitialization());

	// Initialize the visualization
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();

	// The User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	//process macro or start UI session
	//
	if (!ui){
	//batch mode.
	G4String command = "/control/execute ";
	G4String filename = argv[1];
	UImanager->ApplyCommand(command+filename);
	} 
	else {
	// interactive mode.
	UImanager->ApplyCommand("/control/execute init_vis.mac");
	ui->SessionStart();
	delete ui;
	}

	delete visManager;
	delete runManager;

}
