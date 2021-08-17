// trianglesActionInitialization.cc
// Implemantation of tha trianglesActionInitialization class.
// The definitions should be found in (source directory of this project)/src/trianglesActionInitialization.hh
//
#include "trianglesActionInitialization.hh"

#include "trianglesPrimaryGeneratorAction.hh"
#include "trianglesRunAction.hh"
#include "trianglesEventAction.hh"

//The constructor
trianglesActionInitialization::trianglesActionInitialization() : G4VUserActionInitialization() { }

//The destructor
trianglesActionInitialization::~trianglesActionInitialization() { }

void trianglesActionInitialization::BuildForMaster() const{
	trianglesEventAction* eventAction = new trianglesEventAction();
	SetUserAction(new trianglesRunAction(eventAction));
}

void trianglesActionInitialization::Build() const{
	// Set the default properties of the particles.
	SetUserAction(new trianglesPrimaryGeneratorAction);

	trianglesEventAction* eventAction = new trianglesEventAction();
	SetUserAction(eventAction);

	SetUserAction(new trianglesRunAction(eventAction));

}
