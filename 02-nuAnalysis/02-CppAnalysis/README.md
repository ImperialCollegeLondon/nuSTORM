# nuSTORM -- nuAnalis

The code in this directory tree provides a c++ analysis skeleton that can be used to read the nuSIM ntuples (see [nuSTORM](https://www.nustorm.org/trac/).  Documentation for the nuAnalysis package will be posted on the nuSTORM wiki.

## To set up and run:
Execute "startup.bash" from this directory (i.e. run the bash command "source startup.bash").  This will:
  * Set the "nuAnalysisPATH" to this directory; and
  * Set "nuAnalysisCPATH" so that g++ can see the 01-Code directory.  The scripts in "02-Tests" may then be run using the compile scriptis included in 02-Tests.  A sample user analysis programme is included in "03-Skeleton"

## Directories:
 * C++ classes and "library" code stored in "01-Code"
 * Test scripts stored in "02-Tests"
 * Integration tests are stored in "03-Integration-Test"
 * Sample data ntuples are stored in "31-Data"

## Dependencies:
 * g++11 and assumes ROOT is installed.

## Docker Support

A Dockerfile is provided for the convenient writing, building, and testing of code without the need to install the entire ROOT project.
The latest image published by the ROOT project on DockerHub is utilised as the base which simplifies and speeds up getting started.
The Dockerfile is at the root of this 02-CppAnal folder and an image can be built simply by running `docker build -t nustorm_docker .` from within this directory.
This image can then be used and explored by running it with: `docker run --rm -it nustorm_docker` which opens a `bash` shell inside the image. 
Note that any changes made within the image are lost once you exit.

The Dockerfile can also be used to run tests. A script to assist with this is provided - `run_tests_in_docker.sh`
This builds a fresh Dockerfile if required and runs the integration test within it.

### Dependencies:
 * A suitable Docker daemon must be installed. See [Docker](https://www.docker.com/) for instructions

## Style Guide:
This project aspires to follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html)

## History
 * 06 December 2021:  First version.
 * 25 November 2022:  Add Docker Support.

