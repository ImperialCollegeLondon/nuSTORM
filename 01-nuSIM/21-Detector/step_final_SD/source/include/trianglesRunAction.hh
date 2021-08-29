// trianglesRunAction.hh
// Definition of the trianglesRunManager class.
// This class defines the actions that must be taken before, during and at the end of a run.
// However, practically speaking, everything is related to the analysis.

#ifndef trianglesRunAction_h
#define trianglesRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class trianglesEventAction;

class G4Run;

class trianglesRunAction : public G4UserRunAction{
	public:
		trianglesRunAction(trianglesEventAction* eventAction);
		virtual ~trianglesRunAction();

		virtual void BeginOfRunAction(const G4Run*);
		virtual void EndOfRunAction(const G4Run*);

	private:
		trianglesEventAction* fEventAction;
};

#endif
