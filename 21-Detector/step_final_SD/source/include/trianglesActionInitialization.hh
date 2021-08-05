// trianglesActionInitialization.hh
// A brief defenition of the trianglesActionInitialization class.
// The implementation should be found in (source directory of this project)/src/trianglesActionInitialization.cc
//
// This class defines(?) the actions an application user can take, and its initializations.
//
#ifndef trianglesActionInitialization_h
#define trianglesActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class trianglesActionInitialization : public G4VUserActionInitialization{
	public:
		trianglesActionInitialization();
		virtual ~trianglesActionInitialization();

		virtual void BuildForMaster() const;
		virtual void Build() const;

};

#endif
