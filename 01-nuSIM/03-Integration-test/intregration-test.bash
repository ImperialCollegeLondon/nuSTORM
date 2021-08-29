#!/bin/bash
#
#..  Integration test script.  Runs each test and creates integrated log
#    file.  Check of log file constitutes the test.  Can not do "did" as
#    some scripts create output with randomly generated values which will
#    change execution to execution.
#
#..  Assumes "startup.bash" has been run.

cd $nuSIMPATH
02-Tests/MuonConstTst.py
02-Tests/MuonDecayTst.py
02-Tests/MuonDecayTrncTst.py
02-Tests/nuSTORMPrdStrghtTst.py
02-Tests/NeutrinoEventInstanceTst.py
02-Tests/MuonBeam4CoolingDemoTst.py
02-Tests/SimulationTst.py
02-Tests/RunSimulation.py
02-Tests/FluxCalcOutline.py
