# nuSTORM

The code in this repository provides a simulation of neutrino production in [nuSTORM](https://www.nustorm.org/trac/).  The contents of the present simulation and the data structure that is produced is summarised in [nuSIM-2021-01](https://www.nustorm.org/trac/raw-attachment/wiki/Software-and-computing/Documentation/2021/nuSIM-doc-01.pdf).

## To set up and run:
A guide to setting up and running the simulation package "nuSIM" is given in [nuSIM-2021-02](https://www.nustorm.org/trac/raw-attachment/wiki/Software-and-computing/Documentation/2021/nuSIM-doc-02.pdf)
Execute "startup.bash" from this directory (i.e. run the bash command "source startup.bash").  This will:
  * Set up "nuSIMPATH"; and
  * Add "01-Code" to the PYTHONPATH.  The scripts in "02-Tests" may then be run with the command "python 02-Tests/<filename>.py".

## Directories:
 * Python classes and "library" code stored in "01-Code"
 * Test scripts stored in "02-Tests"
 * Integration tests are stored in "03-Integration-Test"
 * Parameters to control the run conditions are stored in "11-Parameters"

Rudimentary, but, goal is one test script per class/package file in 01-Code.

## Dependencies:
 * Code and test scripts assume Python 3.  
 * Test scripts assume code directory (01-Code) is in PYTHON path.  A bash script "startup.bash" is provided to update the PYTHON path.
