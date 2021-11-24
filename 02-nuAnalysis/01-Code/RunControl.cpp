/*  Run control class implementation file */

#include <cstdlib>
#include <iostream>
#include <filesystem>

#include "RunControl.hpp"

RunControl* RunControl::instance = NULL;

RunControl* RunControl::getInstance() {
    if (instance == NULL) {
      
        instance     = new RunControl();
    }
    return instance;
}
  
RunControl::RunControl(bool Dbg,
		       std::string Rfn,
		       std::string Cdn){
  Debug        = Dbg;
  ROOTfilename = Rfn;
  CHAINdirname = Cdn;
}

void RunControl::print() {
  std::cout << "RunControl initialised with parameters:" << std::endl;
  std::cout << "                Debug: " << Debug        << std::endl;
  std::cout << "       ROOT file name: " << ROOTfilename << std::endl;
  std::cout << " ROOT chain directory: " << CHAINdirname << std::endl;
}

void DumpHelp() {
  std::cout << "RunControl help: explanation of arguments:" << std::endl;
  std::cout << "\t -h \t\t Generates this help printout to show use of flags etc." << std::endl;
  std::cout << "\t -d \t\t Sets debug flag: RunControl::Debug = true" << std::endl;
  std::cout << "\t -f <filename> \t Sets ROOT <filename> to be read (single file)." << std::endl;
  std::cout << "\t\t\t If -f is specified, <filename> must exist or execution is" << std::endl;
  std::cout << "\t\t\t terminated" << std::endl;
  std::cout << "\t -c <dir name> \t Sets directory containing ROOT files to be chained."
	    << std::endl;
  std::cout << "\t\t\t If -c is specified, <dir name> must exist or execution is" << std::endl;
  std::cout << "\t\t\t terminated" << std::endl;
  std::cout << "\t Note that -f and -c can be used together, all files will be read."
	    << std::endl;

  std::exit(EXIT_SUCCESS);
}

void RunControl::ParseArgs(int nArgs, char *ArgV[]){
  char *Arg;

  Debug     = false;
  FileFlag  = false;
  ChainFlag = false;
  for (int i = 0 ; i < nArgs ; i++) {
    Arg = ArgV[i];
    if ( std::strcmp(Arg, "-d") == 0 )
      Debug = true;
    else if ( std::strcmp(Arg, "-f") == 0 ) {
      FileFlag     = true;
      ROOTfilename = ArgV[i+1]; }
    else if ( std::strcmp(Arg, "-c") == 0 ) {
      ChainFlag    = true;
      CHAINdirname = ArgV[i+1]; }
    else if ( std::strcmp(Arg, "-h") == 0 ) 
      DumpHelp();
  }
  if ( Debug == true) {
    std::cout << "Debug: parsed " << nArgs << " arguments \n" ;
    for (int i = 0 ; i < nArgs ; i++) {
      std::cout << "    ----> Argument: " << i 
		<< " value: " << ArgV[i] << "\n" ;
    }
  }
  bool File = std::filesystem::exists(ROOTfilename);
  bool Dir  = std::filesystem::exists(CHAINdirname);
  if ( FileFlag and !File ) {
    std::cout << " Error! File "
	      << ROOTfilename << " does not exist; stop." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ( ChainFlag and !Dir ) {
    std::cout << " Error! Directory "
	      << CHAINdirname << " does not exist; stop." << std::endl;
    std::exit(EXIT_FAILURE);
  }
}
