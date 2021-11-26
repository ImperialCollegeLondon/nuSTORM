/* Tests of RunControl class */

#include <iostream>

#include "RunControl.hpp"

int main(int nArgs, char *ArgV[]){

  std::cout << "Start tests of RunControl class:" << std::endl;

  int i = 0;
  // Test of instantiation of singleton class:
  std::cout << "    Test " << i << ": test singleton nature of class:"
	    << std::endl;

  // Initialise run control singleton class:
  RunControl* RC = RunControl::getInstance();
  std::cout << "    ----> Initial RunControl instance: " << RC << std::endl;
  RC->print();
  RunControl* RC1 = RunControl::getInstance();
  std::cout << "    ----> Second RunControl instance: " << RC1 << std::endl;
  if (RC != RC1) {
    std::cout << "    ----> FAILED!  Not a singleton!" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  else {
    std::cout << "        ----> SUCCESS!  A singleton!" << std::endl;
  }
  i++;
  // Test parsing of input arguments:
  std::cout << "    Test " << i << ": test parsing of input arguments:"
	    << std::endl;
  std::cout << "    ----> Number of input arguments: "
	    << nArgs << std::endl;

  std::string Hlp = "-h";
  const char* cHlp = Hlp.c_str();
  char *ArgVPlus[2];
  int nArgsPlus;
  if ( nArgs == 1 ){
    std::cout << "    ----> No arguments given, set -h, help flag."
	      << std::endl;
    ArgVPlus[0] = ArgV[0]; 
    ArgVPlus[1] = const_cast<char *>(cHlp);
    std::cout << "; ArgPlus: " << ArgVPlus[nArgs] << std::endl;
    nArgsPlus = nArgs + 1;
    std::cout << "        ----> Updated: nArgs, ArgV: "
	      << nArgsPlus << ": " << ArgVPlus << std::endl;
    RC->ParseArgs(nArgsPlus, ArgVPlus);
  }
  else {
    RC->ParseArgs(nArgs, ArgV);
  }
  std::cout << "    ----> RunControl parameters:" << std::endl;
  RC->print();
}
